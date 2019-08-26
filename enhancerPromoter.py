#!/usr/bin/python

'''
The MIT License (MIT)

Copyright (c) 2016 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''

#pythonTemplate.py <- change to title of your script
#130801 <- date
#Name 


#Description:

#This is a generic python template that has functions from utils.py imported and can be used on CFCE1



#================================================================================
#=============================DEPENDENCIES=======================================
#================================================================================



import sys
import os
print "Using python version %s" % sys.version


#importing utils package
#add locations of files and global parameters in this section
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = whereAmI + '/'

sys.path.append(pipeline_dir)

import argparse
import cPickle
import utils
import pipeline_dfci
import subprocess

import string
import tempfile
import zlib
import numpy
import re
import time
from distutils.spawn import find_executable
from collections import defaultdict


#================================================================================
#============================GLOBAL PARAMETERS===================================
#================================================================================

#add locations of files and global parameters in this section
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = whereAmI + '/'


bamliquidator_path = 'bamliquidator_batch.py'

#using a paramater dictionary in liue of a yaml or json for now



paramDict = {'cpgPath': '/storage/cylin/grail/projects/mycn_resub/mycn/beds/hg19_cpg_islands.bed',

             
             }



#================================================================================
#===================================CLASSES======================================
#================================================================================

#user defined classes here

#================================================================================
#=================================FUNCTIONS======================================
#================================================================================


def loadAnnotFile(genome,tss_window,geneList=[],skip_cache=False):
    """
    load in the annotation and create a startDict and tss collection for a set of refseq IDs a given genome
    """
    genomeDict = {
        'HG18': 'annotation/hg18_refseq.ucsc',
        'MM9': 'annotation/mm9_refseq.ucsc',
        'MM10': 'annotation/mm10_refseq.ucsc',
        'HG19': 'annotation/hg19_refseq.ucsc',
        'HG19_RIBO': 'annotation/hg19_refseq.ucsc',
        'RN4': 'annotation/rn4_refseq.ucsc',
        'RN6': 'annotation/rn6_refseq.ucsc',
        'HG38': 'annotation/hg38_refseq.ucsc',
        }

    genomeDirectoryDict = {
        'HG19':'/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/',
        'RN6':'/storage/cylin/grail/genomes/Rattus_norvegicus/UCSC/rn6/Sequence/Chromosomes/',
        'MM9':'/storage/cylin/grail/genomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/',
        'MM10':'/storage/cylin/grail/genomes/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/',
        'HG38': '/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/',
        }

    mouse_convert_file = '%s/annotation/HMD_HumanPhenotype.rpt' % (whereAmI)

    #making a dictionary for mouse to human conversion
    mouse_convert_dict = defaultdict(str)
    
    mouse_convert_table = utils.parseTable(mouse_convert_file,'\t')
    for line in mouse_convert_table:
        mouse_convert_dict[line[4]] = line[0]
            
    genomeDirectory = genomeDirectoryDict[string.upper(genome)]


    #making a chrom_dict that is a list of all chroms with sequence
    chrom_list = utils.uniquify([name.split('.')[0] for name in os.listdir(genomeDirectory) if len(name) >0])
    
    annotFile = whereAmI + '/' + genomeDict[string.upper(genome)]


    if not skip_cache:
        # Try loading from a cache, if the crc32 matches
        annotPathHash = zlib.crc32(annotFile) & 0xFFFFFFFF  # hash the entire location of this script
        annotFileHash = zlib.crc32(open(annotFile, "rb").read()) & 0xFFFFFFFF

        cache_file_name = "%s.%s.%s.cache" % (genome, annotPathHash, annotFileHash)

        cache_file_path = '%s/%s' % (tempfile.gettempdir(), cache_file_name)

        if os.path.isfile(cache_file_path):
            # Cache exists! Load it!
            try:
                print('\tLoading genome data from cache.')
                with open(cache_file_path, 'rb') as cache_fh:
                    cached_data = cPickle.load(cache_fh)
                    print('\tCache loaded.')
                return cached_data
            except (IOError, cPickle.UnpicklingError):
                # Pickle corrupt? Let's get rid of it.
                print('\tWARNING: Cache corrupt or unreadable. Ignoring.')
        else:
            print('\tNo cache exists: Loading annotation (slow).')


    # We're still here, so either caching was disabled, or the cache doesn't exist

    startDict = utils.makeStartDict(annotFile, geneList)
    tssLoci =[]
    if geneList==[]:
        geneList = startDict.keys()
    for gene in geneList:
        tssLoci.append(utils.makeTSSLocus(gene,startDict,tss_window,tss_window))

    tssCollection = utils.LocusCollection(tssLoci,50)

    if not skip_cache:
        print('Writing cache for the first time.')
        with open(cache_file_path, 'wb') as cache_fh:
            cPickle.dump((startDict, tssCollection), cache_fh, cPickle.HIGHEST_PROTOCOL)

    return startDict, tssCollection, genomeDirectory, chrom_list, mouse_convert_dict




def splitRegions(inputGFF,tssCollection):


    #if even a single coordinate is shared with the +/-1kb 
    splitGFF = []
    debugCount = 0
    for line in inputGFF:

        chrom = line[0]
        regionID = line[1]
        lineLocus = utils.Locus(line[0],line[3],line[4],'.')

        overlappingLoci = tssCollection.getOverlap(lineLocus)
        if len(overlappingLoci) > 0: #case where a tss Overlap
            #identify the parts of the line locus that are contained
            localTSSCollection = utils.LocusCollection(overlappingLoci,50)
            overlappingCoords = lineLocus.coords()
            for tssLocus in overlappingLoci:
                overlappingCoords += tssLocus.coords()

            overlappingCoords = utils.uniquify(overlappingCoords)
            overlappingCoords.sort()

            #you need to hack and slash add 1 to the last coordinate of the overlappingCoords
            overlappingCoords[-1] +=1

            i = 0
            regionTicker = 1
            while i < (len(overlappingCoords)-1):
                start = int(overlappingCoords[i])
                stop = int(overlappingCoords[(i+1)])-1
                if (stop - start) < 50: #this eliminates really tiny regions
                    i+=1
                    continue
                splitLocus = utils.Locus(chrom,start+1,stop,'.')


                if lineLocus.overlaps(splitLocus): #has to be a mycn site
                    newID = '%s_%s' % (regionID,regionTicker)
                    tssStatus = 0
                    if localTSSCollection.getOverlap(splitLocus):
                        tssStatus = 1
                    splitGFFLine = [chrom,newID,newID,start,stop,'','.',tssStatus,newID]

                    splitGFF.append(splitGFFLine)
                    regionTicker+=1
                i+=1
        else:
            line[7] = 0
            splitGFF.append(line)




    return splitGFF




def mapBams(bamFileList,splitGFFPath,analysisName,mappedFolder):

    print("MAPPING TO THE FOLLOWING BAMS:")

    for bamFile in bamFileList:
        print(bamFile)
        bamFileName = bamFile.split('/')[-1]

        # MAPPING TO THE STITCHED GFF
        mappedOut1Folder = '%s%s_%s_MAPPED' % (mappedFolder, analysisName, bamFileName)
        mappedOut1File = '%s%s_%s_MAPPED/matrix.txt' % (mappedFolder, analysisName, bamFileName)
        if utils.checkOutput(mappedOut1File, 0.2, 0.2):
            print("FOUND %s MAPPING DATA FOR BAM: %s" % (splitGFFPath, mappedOut1File))
        else:
            cmd1 = bamliquidator_path + " --sense . -e 200 --match_bamToGFF -r %s -o %s %s" % (splitGFFPath, mappedOut1Folder, bamFile)
            print(cmd1)

            os.system(cmd1)
            if utils.checkOutput(mappedOut1File,0.2,5):
                print("SUCCESSFULLY MAPPED TO %s FROM BAM: %s" % (splitGFFPath, bamFileName))
            else:
                print("ERROR: FAILED TO MAP %s FROM BAM: %s" % (splitGFFPath, bamFileName))
                sys.exit()

    print('BAM MAPPING COMPLETED NOW MAPPING DATA TO REGIONS')
    
    #now we make a signal table
    #set up the table using the first bam
    if len(bamFileList) > 1:

        #set up the first pass at the table
        signalTable = [['REGION_ID','locusLine'] + [name.split('/')[-1] for name in bamFileList]]
        bamFileName = bamFileList[0].split('/')[-1]
        mappedTable = utils.parseTable( '%s%s_%s_MAPPED/matrix.txt' % (mappedFolder, analysisName, bamFileName),'\t')
        for i in range(1,len(mappedTable)):
            signalTable.append(mappedTable[i])

        for bamFile in bamFileList[1:]:
            bamFileName = bamFile.split('/')[-1]
            
            mappedTable = utils.parseTable( '%s%s_%s_MAPPED/matrix.txt' % (mappedFolder, analysisName, bamFileName),'\t')
            
            for i in range(1,len(mappedTable)):
                mapSignal = mappedTable[i][2]
                signalTable[i].append(mapSignal)
    else:
        bamFileName = bamFileList[0].split('/')[-1]
        signalTable = utils.parseTable( '%s%s_%s_MAPPED/matrix.txt' % (mappedFolder, analysisName, bamFileName),'\t')

    return(signalTable)



def makeAverageTable(outputFolder,analysisName,useBackground = False):
    '''
    makes a signal table that is the average background subtracted signal for each region
    if background is present, will zero out regions before trying to take average. i.e. no negative regions allowed
    '''
    
    #first the easy case with no background
    if not useBackground:
        signalTablePath = '%s%s_signal_table.txt' % (outputFolder,analysisName)
        signalTable = utils.parseTable(signalTablePath,'\t')

        averageTable = [['GENE_ID','locusLine','%s_signal' % (analysisName)]]

        for line in signalTable[1:]:
            newLine = line[0:2]
            avgSignal = round(numpy.mean([float(x) for x in line[2:]]),4)
            newLine.append(avgSignal)
            averageTable.append(newLine)
    #now the condition w/ background
    else:
        signalTablePath = '%s%s_signal_table.txt' % (outputFolder,analysisName)
        signalTable = utils.parseTable(signalTablePath,'\t')

        controlTablePath = '%s%s_control_signal_table.txt' % (outputFolder,analysisName)
        controlTable = utils.parseTable(controlTablePath,'\t')

        averageTable = [['GENE_ID','locusLine','%s_signal' % (analysisName)]]

        #checking to make sure the # of backgrounds = number of signal bams
        #otherwise throw an error
        signal_n_col = len(signalTable[0])
        control_n_col = len(controlTable[0])

        if signal_n_col != control_n_col:
            print('ERROR: MUST PROVIDE SAME NUMBER OF CONTROL BAMS')
            sys.exit()

        signal_n_rows = len(signalTable)
        control_n_rows = len(controlTable)

        if signal_n_rows != control_n_rows:
            print('ERROR: MAPPED FILES ARE NOT THE SAME LENGTH')
            sys.exit()

        for i in range(1,len(signalTable)):
            signalLine = signalTable[i]
            controlLine = controlTable[i]
            if signalLine[0:2] != controlLine[0:2]:
                print('ERROR: REGIONS ON LINE %s DO NOT CORRESPOND' % (i))
                sys.exit()

            newLine = signalLine[0:2]

            signal_values = [float(x) for x in signalLine[2:]]
            control_values = [float(x) for x in controlLine[2:]]

            subtracted_values = [signal_values[x] - control_values[x] for x in range(len(signal_values))]
            subtracted_values = [max(0,x) for x in subtracted_values] # now make negative numbers 0
            avgSignal = round(numpy.mean(subtracted_values),4)
            newLine.append(avgSignal)
            averageTable.append(newLine)

    return averageTable



def makePeakTable(paramDict,splitGFFPath,averageTablePath,startDict,geneList,genomeDirectory,tss_window,distal_window,tads_path=''):
    
    '''
    makes the final peak table with ebox info
    '''

    peakTable = [['REGION_ID','CHROM','START','STOP','LENGTH','TSS','CPG','CPG_FRACTION','GC_FREQ','SIGNAL','CANON_EBOX_COUNT','NON_CANON_EBOX_COUNT','TOTAL_EBOX_COUNT','OVERLAPPING_GENES','PROXIMAL_GENES']]
    
    print('LOADING PEAK REGIONS')
    peakGFF = utils.parseTable(splitGFFPath,'\t')

             
    print('LOADING BINDING DATA')
    signalTable = utils.parseTable(averageTablePath,'\t')

    print('LOADING CPGS ISLANDS')
    cpgBed = utils.parseTable(paramDict['cpgPath'],'\t')
    cpgLoci = []
    for line in cpgBed:
        cpgLoci.append(utils.Locus(line[0],line[1],line[2],'.',line[-1]))
    cpgCollection = utils.LocusCollection(cpgLoci,50)

    print("MAKING TSS COLLECTIONS")
    if len(geneList) == 0:
        geneList = startDict.keys()

    tss_prox_loci = []
    tss_distal_loci = []
    for refID in geneList:
        tss_prox_loci.append(utils.makeTSSLocus(refID,startDict,tss_window,tss_window))
        tss_distal_loci.append(utils.makeTSSLocus(refID,startDict,distal_window,distal_window))

    #make a 1kb flanking and 50kb flanking collection
    tss_prox_collection = utils.LocusCollection(tss_prox_loci,50)
    tss_distal_collection = utils.LocusCollection(tss_distal_loci,50)

    if len(tads_path) > 0:
        print('LOADING TADS FROM %s' % (tads_path))
        tad_collection = utils.importBoundRegion(tads_path,'tad')
        use_tads = True
        
        #building a tad dict keyed by tad ID w/ genes in that tad provided
        tad_dict = defaultdict(list)
        for tss_locus in tss_prox_loci:
            overlapping_tads = tad_collection.getOverlap(tss_locus,'both')
            for tad_locus in overlapping_tads:
                tad_dict[tad_locus.ID()].append(tss_locus.ID())

    else:
        use_tads = False

    print('CLASSIFYING PEAKS')
    ticker = 0

    no_tad_count = 0
    for i in range(len(peakGFF)):
        if ticker%1000 == 0:
            print(ticker)
        ticker +=1

        #getting the particulars of the region
        gffLine = peakGFF[i]
        peakID = gffLine[1]
        chrom = gffLine[0]
        start = int(gffLine[3])
        stop = int(gffLine[4])
        lineLocus = utils.Locus(chrom,start,stop,'.',peakID)

        #getting the mapped signal
        signalLine = signalTable[(i+1)]        
        signalVector = [float(x) for x in signalLine[2:]]
        


        #setting up the new line
        newLine = [peakID,chrom,start,stop,lineLocus.len()]

        #get the tss status from the gff itself (we are able to do this nicely from the split gff code earlier
        newLine.append(gffLine[7])
        
        #check cpg status
        if cpgCollection.getOverlap(lineLocus,'both'):
            newLine.append(1)
        else:
            newLine.append(0)

        #now do fractional cpgOverlap
        overlappingCpGLoci = cpgCollection.getOverlap(lineLocus,'both')
        overlappingBases = 0
        for locus in overlappingCpGLoci:
            cpgStart = max(locus.start(),lineLocus.start())
            cpgEnd = min(locus.end(),lineLocus.end())
            overlappingBases += (cpgEnd-cpgStart)
        overlapFraction = float(overlappingBases)/lineLocus.len()
        
        newLine.append(round(overlapFraction,2))

        #now get the seq
        lineSeq = string.upper(utils.fetchSeq(genomeDirectory,chrom,start,stop,True))
        if len(lineSeq) == 0:
            print('UH OH')
            print(lineSeq)
            print(gffLine)
            print(i)
            print(chrom)
            print(start)
            print(stop)
            sys.exit()
        
        gcFreq = float(lineSeq.count('GC') + lineSeq.count('CG'))/len(lineSeq)
        newLine.append(gcFreq)
            
        #this is where we add the ChIP-Seq signal
        newLine += signalVector


        eboxMatchList = re.findall('CA..TG',lineSeq)
        if len(eboxMatchList) == 0:
            newLine += [0]*3
        else:
            totalCount = len(eboxMatchList)
            canonCount = eboxMatchList.count('CACGTG')
            otherCount = totalCount - canonCount
            newLine += [canonCount,otherCount,totalCount]

        #now find the overlapping and proximal genes
        #here each overlapping gene the tss prox locus overlaps the peak

        if use_tads:

            tad_loci = tad_collection.getOverlap(lineLocus,'both')

            tad_id_list = [tad_locus.ID() for tad_locus in tad_loci]
            tad_genes = []
            for tad_id in tad_id_list:
                tad_genes+=tad_dict[tad_id]
            if len(tad_genes) == 0:
                #print('no tad for this region')
                #print(gffLine)
                no_tad_count+=1
        else:
            tad_genes=[]



        if len(tad_genes) >0:
            overlappingGenes = [startDict[locus.ID()]['name'] for locus in tss_prox_collection.getOverlap(lineLocus,'both') if tad_genes.count(locus.ID()) > 0]
            proximalGenes = [startDict[locus.ID()]['name'] for locus in tss_distal_collection.getOverlap(lineLocus,'both') if tad_genes.count(locus.ID()) > 0]        
            # print('linked peak to tad genes')
            # print([startDict[x]['name'] for x in tad_genes])
            # print(tad_id_list)
            # print(gffLine)
            # print(overlappingGenes)
            # print(proximalGenes)
        else:
            overlappingGenes = [startDict[locus.ID()]['name'] for locus in tss_prox_collection.getOverlap(lineLocus,'both')]
            proximalGenes = [startDict[locus.ID()]['name'] for locus in tss_distal_collection.getOverlap(lineLocus,'both')]        

        overlappingGenes = utils.uniquify(overlappingGenes) 
        #here the tss 50kb locus overlaps the peak 
        #overlap takes priority over proximal
        proximalGenes = [gene for gene in proximalGenes if overlappingGenes.count(gene) == 0]
        proximalGenes = utils.uniquify(proximalGenes)


        overlappingString = string.join(overlappingGenes,',')
        proximalString = string.join(proximalGenes,',')
        

        newLine += [overlappingString,proximalString]


        peakTable.append(newLine)

    print('Out of %s regions, %s were assigned to at least 1 tad' % (len(peakTable),no_tad_count))
    return peakTable


def makeGeneTable(peakTable,analysisName):

    '''
    takes the peak table and makes a gene centric table

    '''


    geneDict = {}

    geneTable = [['GENE','%s_TSS_SIGNAL' % (analysisName),'%s_DISTAL_SIGNAL' % (analysisName)]]

    #now iterate through the table
    for line in peakTable[1:]:
        regionLength = int(line[4])

        signal = float(line[9]) * regionLength


        #genes where this particular peak overlaps the tss prox window
        #where there are both overlap and proximal meet
        if len(line) == 15:
            overlapGeneList = [gene for gene in line[-2].split(',') if len(gene) > 0]
            if overlapGeneList.count('107'):
                print(line)
                sys.exit()
            for overlapGene in overlapGeneList:
                if geneDict.has_key(overlapGene) == False:
                    geneDict[overlapGene] = {'tss':0.0,'distal':0.0}
                #there can be a nasty 1 overlap case where the region might overlap by the overlapping gene list, but not be real
                if int(line[5]) == 1:
                    geneDict[overlapGene]['tss'] += signal
                else: #this is the case where the mycn site is just outside of the promoter or overlapping the gene locus/body these are rar
                    geneDict[overlapGene]['distal'] += signal

            proximalGeneList = [gene for gene in line[-1].split(',') if len(gene) > 0]
            for proximalGene in proximalGeneList:
                if geneDict.has_key(proximalGene) == False:
                    geneDict[proximalGene] = {'tss':0.0,'distal':0.0}
                if int(line[5]) == 0:
                    geneDict[proximalGene]['distal'] += signal
        #where there's just overlap
        if len(line) == 14:
            overlapGeneList = [gene for gene in line[-1].split(',') if len(gene) > 0]
            if overlapGeneList.count('107'):
                print(line)
                sys.exit()
            for overlapGene in overlapGeneList:
                if geneDict.has_key(overlapGene) == False:
                    geneDict[overlapGene] = {'tss':0.0,'distal':0.0}
                #there can be a nasty 1 overlap case where the region might overlap by the overlapping gene list, but not be real
                if int(line[5]) == 1:
                    geneDict[overlapGene]['tss'] += signal
                else: #this is the case where the mycn site is just outside of the promoter or overlapping the gene locus/body these are rar
                    geneDict[overlapGene]['distal'] += signal



    geneList = geneDict.keys()
    geneList = utils.uniquify(geneList)
    geneList.sort()

    for gene in geneList:
        newLine = [gene]
        newLine.append(geneDict[gene]['tss'])
        newLine.append(geneDict[gene]['distal'])
        geneTable.append(newLine)

    return geneTable


def callRWaterfall(geneTablePath,outputFolder,analysisName,top):
    
    '''
    function to call the Rscript and to wait until the .cls and .gct files are created
    returns the paths
    '''

    rBashFilePath = '%s%s_R_plotting.sh' % (outputFolder,analysisName)
    rBashFile = open(rBashFilePath,'w')
    rBashFile.write('#!/usr/bin/bash\n\n')

    rCmd = 'Rscript %senhancerPromoter_waterfall.R %s %s %s %s' % (pipeline_dir,geneTablePath,outputFolder,analysisName,top)
    rBashFile.write(rCmd)
    rBashFile.close()
    print('writing R plotting command to disk and calling %s' %(rBashFilePath))
    os.system('bash %s' % (rBashFilePath))

    #now check for the .cls output
    clsPath = '%s%s_top_all.cls' % (outputFolder,analysisName)

    if utils.checkOutput(clsPath,0.5,5):
        return 
    else:
        print('ERROR: UNABLE TO SUCCESFULLY DETECT R SCRIPT OUTPUT AT %s' % (clsPath))
        sys.exit()


def callGSEA(outputFolder,analysisName,top,analysis_type ='enhancer_vs_promoter',use_top=True):

    '''
    runs C2 GSEA
    '''

    #figure out the suffix for gct and cls files
    analysis_dict = {'enhancer_vs_promoter':['','#PROMOTER_versus_DISTAL'],
                   'total_contribution':['_total_contrib','#SIGNAL_versus_BACKGROUND'],
                   }
    
    if analysis_dict.has_key(analysis_type) == False:
        print('Error: please use one of the following supported analysis types')
        print(analysis_dict.keys())
        sys.exit()

    suffix = analysis_dict[analysis_type][0]

    gseaPath = '/storage/cylin/home/cl6/gsea2-3.0_beta_2.jar'
    gmxPath = '/storage/cylin/grail/annotations/gsea/c2.all.v5.1.symbols.gmt' #C2 set


    gseaBashFilePath = '%s%s_GSEA%s_cmd.sh' % (outputFolder,analysisName,suffix)
    gseaBashFile = open(gseaBashFilePath,'w')

    gseaBashFile.write('#!/usr/bin/bash\n\n')

    gseaBashFile.write('#COMMAND LINE GSEA CALLS FOR %s USING %s COMPARISON\n\n' % (analysisName,string.upper(analysis_type)))
    
    
    #for all
    gctPath = '%s%s_top_all%s.gct' % (outputFolder,analysisName,suffix)
    clsPath = '%s%s_top_all%s.cls' % (outputFolder,analysisName,suffix)
    gseaOutputFolder = utils.formatFolder('%sgsea_top_all_c2%s' % (outputFolder,suffix),True)
    rptLabel = '%s_top_all%s' % (analysisName,suffix)

    gseaCmd_all = 'java -Xmx4000m -cp %s xtools.gsea.Gsea -res %s -cls %s%s -gmx %s -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label %s -metric Diff_of_Classes -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out %s -gui false' % (gseaPath,gctPath,clsPath,analysis_dict[analysis_type][1],gmxPath,rptLabel,gseaOutputFolder)

    gseaBashFile.write(gseaCmd_all)
    gseaBashFile.write('\n')

    if use_top:
        #for top N
        gctPath = '%s%s_top_%s%s.gct' % (outputFolder,analysisName,top,suffix)
        clsPath = '%s%s_top_%s%s.cls' % (outputFolder,analysisName,top,suffix)
        gseaOutputFolder = utils.formatFolder('%sgsea_top_%s_c2%s' % (outputFolder,top,suffix),True)
        rptLabel = '%s_top_%s%s' % (analysisName,top,suffix)

        gseaCmd_top = 'java -Xmx4000m -cp %s xtools.gsea.Gsea -res %s -cls %s%s -gmx %s -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label %s -metric Diff_of_Classes -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out %s -gui false' % (gseaPath,gctPath,clsPath,analysis_dict[analysis_type][1],gmxPath,rptLabel,gseaOutputFolder)

        gseaBashFile.write(gseaCmd_top)
        gseaBashFile.write('\n')

    gseaBashFile.close()
    os.system('bash %s' % (gseaBashFilePath))
    

def detectGSEAOutput(analysisName,outputFolder,top,analysis_type ='enhancer_vs_promoter'):

    '''
    tries to detect the .xls files that show up when GSEA is done running
    '''

    #figure out the suffix for gct and cls files
    analysis_dict = {'enhancer_vs_promoter':['','#PROMOTER_versus_DISTAL','PROMOTER','DISTAL'],
                     'total_contribution':['_total_contrib','#SIGNAL_versus_BACKGROUND','SIGNAL','BACKGROUND'],
                   }
    
    if analysis_dict.has_key(analysis_type) == False:
        print('Error: please use one of the following supported analysis types')
        print(analysis_dict.keys())
        sys.exit()

    suffix = analysis_dict[analysis_type][0]

    #first figure out the friggin output folder 
    gseaParentFolder = '%sgsea_top_%s_c2%s/' % (outputFolder,top,suffix)

    for i in range(30):
        folderList = os.listdir(gseaParentFolder)
        #print(folderList)
        
        candidateFolderList = [folder for folder in folderList if folder.count('%s_top_%s%s.Gsea' % (analysisName,top,suffix)) == 1]
        if len(candidateFolderList) > 1:
            print('ERROR: MULTIPLE GSEA OUTPUT FOLDERS DETECTED FOR %s WITH TOP %s GENES' % (analysisName,string.upper(str(top))))
            sys.exit()
        elif len(candidateFolderList) == 0:
            time.sleep(10)
        elif len(candidateFolderList) == 1:
            candidateFolder = '%sgsea_top_%s_c2%s/%s/' % (outputFolder,top,suffix,candidateFolderList[0])

    print('USING %s AS CANDIDATE GSEA FOLDER' % (candidateFolder))
    timeStamp = candidateFolder.split('.')[-1][:-1]
    print(timeStamp)
    #now that you have the candidate folder find the friggen xls files
    
    class_1 = analysis_dict[analysis_type][2]
    class_2 = analysis_dict[analysis_type][3]
    #for promoter
    class1TablePath = '%sgsea_report_for_%s_%s.xls' % (candidateFolder,class_1,timeStamp)
    class2TablePath = '%sgsea_report_for_%s_%s.xls' % (candidateFolder,class_2,timeStamp)
    print(class1TablePath)
    print(class2TablePath)
    #now check em
    if utils.checkOutput(class1TablePath,0.5,30):
        print('FOUND %s OUTPUT AT %s' % (class_1,class1TablePath))
        if utils.checkOutput(class2TablePath,0.5,30):
            print('FOUND %s OUTPUT AT %s' % (class_2,class2TablePath))
            return class1TablePath,class2TablePath
    else:
        print('ERROR: UNABLE TO FIND GSEA OUTPUT')

def callR_GSEA(class1TablePath,class2TablePath,outputFolder,analysisName,top):
    
    '''
    function to call the Rscript and to wait until the .cls and .gct files are created
    returns the paths
    '''

    rBashFilePath = '%s%s_R_gsea.sh' % (outputFolder,analysisName)
    rBashFile = open(rBashFilePath,'w')
    rBashFile.write('#!/usr/bin/bash\n\n')

    rCmd = 'Rscript %senhancerPromoter_gsea.R %s %s %s %s %s' % (pipeline_dir,class1TablePath,class2TablePath,outputFolder,analysisName,top)
    rBashFile.write(rCmd)
    rBashFile.close()
    print('writing R plotting command to disk and calling %s' %(rBashFilePath))
    os.system('bash %s' % (rBashFilePath))

    #now check for the nes output
    nesPath = '%s%s_top_%s_nes.txt' % (outputFolder,analysisName,top)

    if utils.checkOutput(nesPath,0.5,5):
        return 
    else:
        print('ERROR: UNABLE TO SUCCESFULLY DETECT R SCRIPT OUTPUT AT %s' % (nesPath))
        sys.exit()

#================================================================================
#===============================MAIN RUN=========================================
#================================================================================


#write the actual script here


def main():

    '''
    main run method for enhancer promoter contribution tool
    '''

    parser = argparse.ArgumentParser(usage='%(prog)s [options]')

    # required flags
    parser.add_argument("-b", "--bam", dest="bam", nargs='*',
                        help="Enter a space separated list of .bam files for the main factor", required=True)
    parser.add_argument("-i", "--input", dest="input", type=str,
                        help="Enter .gff or .bed file of regions to analyze", required=True)
    parser.add_argument("-g", "--genome", dest="genome", type=str,
                        help="specify a genome, HG18,HG19,HG38,MM8,MM9,MM10,RN6 are currently supported", required=True)
    

    # output flag
    parser.add_argument("-o", "--output", dest="output", type=str,
                        help="Enter the output folder.", required=True)


    # additional options flags and optional arguments
    parser.add_argument("-a", "--activity", dest="activity", type=str,
                        help="specify a table where first column represents a list of active refseq genes", required=False)

    parser.add_argument("-c", "--control", dest="control", nargs='*',
                        help="Enter a space separated list of .bam files for background. If flagged, will perform background subtraction", required=False)
    parser.add_argument("-t", "--tss", dest="tss",type=int,
                        help="Define the TSS area +/- the TSS. Default is 1kb", required=False, default=1000)
    parser.add_argument("-d", "--distal", dest="distal",type=int,
                        help="Enter a window to assign distal enhancer signal. Default is 50kb", required=False, default=50000)





    parser.add_argument("--other-bams", dest="other", nargs='*',
                        help="enter a space separated list of other bams to map to", required=False)

    parser.add_argument("--name", dest="name", type=str,
                        help="enter a root name for the analysis, otherwise will try to find the name from the input file", required=False)

    parser.add_argument("--top", dest="top", type=int,
                        help="Run the analysis on the top N genes by total signal. Default is 5000", required=False,default=5000)
    parser.add_argument("--tads", dest="tads", type=str,
                        help="Include a .bed of tad regions to restrict enhancer/gene association", required=False,default=None)





    args = parser.parse_args()

    print(args)

    #minimum arguments needed to proceed
    if args.bam and args.input and args.genome and args.output:

        #=====================================================================================
        #===============================I. PARSING ARGUMENTS==================================
        #=====================================================================================

        print('\n\n#======================================\n#===========I. DATA SUMMARY============\n#======================================\n')

        #top analysis subset
        top = args.top


        #input genome
        genome = args.genome.upper()
        print('PERFORMING ANALYSIS ON %s GENOME BUILD' % (genome))
        
        #set of bams
        bamFileList = args.bam

        #bring in the input path
        inputPath = args.input

        #try to get the input name or use the name argument
        if args.name:
            analysisName = args.name
        else:
            analysisName = inputPath.split('/')[-1].split('.')[0]

        print('USING %s AS ANALYSIS NAME' % (analysisName))
        #setting up the output folder
        parentFolder = utils.formatFolder(args.output,True)
        outputFolder = utils.formatFolder('%s%s' % (parentFolder,analysisName),True)

        print('WRITING OUTPUT TO %s' % (outputFolder))


        if inputPath.split('.')[-1] == 'bed':
            #type is bed
            print('input in bed format, converting to gff')
            inputGFF = utils.bedToGFF(inputPath)
        else:
            inputGFF = utils.parseTable(inputPath,'\t')

        
        #the tss window for proximal signal assignment
        tss_window = int(args.tss)

        #the distal window for assigning nearby enhancer signal
        distal_window = int(args.distal)

        #activity path
        if args.activity:
            activityPath = args.activity
            activityTable = utils.parseTable(activityPath,'\t')
            ref_col = 0
            #try to find the column for refseq id
            for i in range(len(activityTable[2])): #use an internal row in case of header
                if str(activityTable[1][i]).count('NM_') > 0 or str(activityTable[1][i]).count('NR_') >0:
                    ref_col = i
            
            #now check for header
            if str(activityTable[0][i]).count('NM_') == 0 and str(activityTable[0][i]).count('NR_') ==0:
                print('REMOVING HEADER FROM GENE TABLE:')
                print(activityTable[0])
                activityTable.pop(0)

            geneList = [line[ref_col] for line in activityTable] # this needs to be REFSEQ NM ID
            print('IDENTIFIED %s ACTIVE GENES' % (len(geneList)))

        else:
            geneList = []

        #check if tads are being invoked
        if args.tads:
            print('LOADING TAD LOCATIONS FROM %s' % (args.tads))
            use_tads = True
            tads_path = args.tads
        else:
            use_tads = False
            tads_path = ''

        print('LOADING ANNOTATION DATA FOR GENOME %s' % (genome))
        
        #important here to define the window
        startDict,tssCollection,genomeDirectory,chrom_list,mouse_convert_dict = loadAnnotFile(genome,tss_window,geneList,True)
        #print(tssCollection.getOverlap(utils.Locus('chr5',171387630,171388066,'.')))
        #sys.exit()

        print('FILTERING THE INPUT GFF FOR GOOD CHROMOSOMES')
        

        print(chrom_list)
        filtered_gff = [line for line in inputGFF if chrom_list.count(line[0]) > 0]
        
        print('%s of INITIAL %s REGIONS ARE IN GOOD CHROMOSOMES' % (len(filtered_gff),len(inputGFF)))

        #=====================================================================================
        #================II. IDENTIFYING TSS PROXIMAL AND DISTAL ELEMENTS=====================
        #=====================================================================================

        print('\n\n#======================================\n#==II. MAPPING TO TSS/DISTAL REGIONS===\n#======================================\n')


        #now we need to split the input region 
        print('SPLITTING THE INPUT GFF USING A WINDOW OF %s' % (tss_window))
        splitGFF = splitRegions(filtered_gff,tssCollection)
        print(len(filtered_gff))
        print(len(splitGFF))

        splitGFFPath = '%s%s_SPLIT.gff' % (outputFolder,analysisName)
        utils.unParseTable(splitGFF,splitGFFPath,'\t')
        print('WRITING TSS SPLIT GFF OUT TO %s' % (splitGFFPath))

        #now you have to map the bams to the gff
        print('MAPPING TO THE SPLIT GFF')
        mappedFolder = utils.formatFolder('%sbam_mapping' % (outputFolder),True)
        
        signalTable = mapBams(bamFileList,splitGFFPath,analysisName,mappedFolder)
        signalTablePath = '%s%s_signal_table.txt' % (outputFolder,analysisName)
        utils.unParseTable(signalTable,signalTablePath,'\t')

        if args.control:
            controlBamFileList = args.control
            controlSignalTable = mapBams(controlBamFileList,splitGFFPath,analysisName,mappedFolder)
            controlSignalTablePath = '%s%s_control_signal_table.txt' % (outputFolder,analysisName)
            utils.unParseTable(controlSignalTable,controlSignalTablePath,'\t')

        #now create the background subtracted summarized average table
        
        print('CREATING AN AVERAGE SIGNAL TABLE')
        averageTable = makeAverageTable(outputFolder,analysisName,useBackground = args.control)
        averageTablePath = '%s%s_average_table.txt' % (outputFolder,analysisName)
        utils.unParseTable(averageTable,averageTablePath,'\t')


        #now load up all of the cpg and other parameters to make the actual peak table

        #first check if this has already been done
        peakTablePath = '%s%s_PEAK_TABLE.txt' % (outputFolder,analysisName)
        if utils.checkOutput(peakTablePath,0.1,0.1):
            print('PEAK TABLE OUTPUT ALREADY EXISTS')
            peakTable = utils.parseTable(peakTablePath,'\t')
        else:
            peakTable = makePeakTable(paramDict,splitGFFPath,averageTablePath,startDict,geneList,genomeDirectory,tss_window,distal_window,tads_path)        
            utils.unParseTable(peakTable,peakTablePath,'\t')

        geneTable = makeGeneTable(peakTable,analysisName)        

        geneTablePath = '%s%s_GENE_TABLE.txt' % (outputFolder,analysisName)
        utils.unParseTable(geneTable,geneTablePath,'\t')

        #if mouse, need to convert genes over
        if genome.count('MM') ==1:
            print('CONVERTING MOUSE NAMES TO HUMAN HOMOLOGS FOR GSEA')
            converted_geneTablePath = '%s%s_GENE_TABLE_CONVERTED.txt' % (outputFolder,analysisName)
        
            converted_geneTable = [geneTable[0]]
            for line in geneTable[1:]:
                converted_name = mouse_convert_dict[line[0]]
                if len(converted_name) >0:
                    converted_geneTable.append([converted_name] + line[1:])

                    utils.unParseTable(converted_geneTable,converted_geneTablePath,'\t')

            geneTablePath = converted_geneTablePath
            geneTable = converted_geneTable

        #=====================================================================================
        #===================================III. PLOTTING ====================================
        #=====================================================================================

        print('\n\n#======================================\n#===III. PLOTTING ENHANCER/PROMOTER===\n#======================================\n')

        #if there are fewer genes in the gene table than the top genes, only run on all
        if len(geneTable)  < int(top):
            print('WARNING: ONLY %s GENES WITH SIGNAL AT EITHER PROMOTERS OR ENHANCERS. NOT ENOUGH TO RUN ANALYSIS ON TOP %s' % (len(geneTable)-1,top))
            top = 0
            use_top =False
        else:
            use_top =True

        #now call the R code
        print('CALLING R PLOTTING SCRIPTS')
        callRWaterfall(geneTablePath,outputFolder,analysisName,top)


        #=====================================================================================
        #==================================IV. RUNNING GSEA===================================
        #=====================================================================================

        print('\n\n#======================================\n#============IV. RUNNING GSEA=========\n#======================================\n')

        #now let's call gsea
        print('RUNNING GSEA ON C2')
        callGSEA(outputFolder,analysisName,top,'enhancer_vs_promoter',use_top)       
        callGSEA(outputFolder,analysisName,top,'total_contribution',use_top)        
        
        if use_top:
            print('DETECTING GSEA OUTPUT FOR TOP %s GENES' % (top))
            #for top by enhancer v promoter metric 
            top_promoterTablePath,top_distalTablePath = detectGSEAOutput(analysisName,outputFolder,top,'enhancer_vs_promoter')
            top_signalTablePath,top_backgroundTablePath = detectGSEAOutput(analysisName,outputFolder,top,'total_contribution')

            print('MAKING NES PLOTS FOR TOP %s GENES' % (top))
            callR_GSEA(top_promoterTablePath,top_distalTablePath,outputFolder,analysisName+'_enhancer_vs_promoter',top)
            callR_GSEA(top_signalTablePath,top_backgroundTablePath,outputFolder,analysisName+'_total_contribution',top)


        print('DETECTING GSEA OUTPUT FOR ALL GENES')
        #for top 
        all_promoterTablePath,all_distalTablePath = detectGSEAOutput(analysisName,outputFolder,'all')

        print('MAKING NES PLOTS FOR ALL GENES')
        callR_GSEA(all_promoterTablePath,all_distalTablePath,outputFolder,analysisName,'all')


        #these files can be parsed to make the NES plot

        #[x for x in fileList if x.count('report_for') == 1and x.count('xls') ==1]
        print('ALL DONE WITH ANALYSIS FOR %s' % (analysisName))
        
main()

