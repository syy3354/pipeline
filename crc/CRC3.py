######################
#
# Core Regulatory Circuits
# Young and Bradner Labs
# Version 1.0
# 140724
#
######################

######################           
# Dependencies       
######################



import os
import sys
import string

whereAmI = os.path.dirname(os.path.realpath(__file__))
pipeline_dir = '%s/' % (string.replace(whereAmI,'crc/','')) # need to set this to where this code is stored
sys.path.append(pipeline_dir)

import utils

import numpy
import scipy
import scipy.stats

import subprocess
import os


from random import randrange
from collections import defaultdict

import networkx as nx
from networkx.algorithms.clique import find_cliques_recursive
import pickle



#================================================================================
#===================================CLASSES======================================
#================================================================================

#user defined classes here

class Genome:
    __chrDict = dict()
    __featureDict = dict()

    #at its core, the genome has attributes of a build name, a fasta directory and an annotation file
    def __init__(self,name,genome_directory,annot_file):
        self._name = name
        self._directory = genome_directory
        self._annot = annot_file

    def name(self):
        return self._name
    
    def directory(self):
        return self._directory

    def annot(self):
        return self._annot

    def addFeature(self,feature,path):

        if feature in self.__featureDict:
            print('WARNING OVERRIDING %s PATH WITH %s' % (feature,path))
        self.__featureDict[feature] = path

    def returnFeature(self,feature):

        #tries to load the selected feature from the feature dictionary

        if feature not in self.__featureDict:
            print('ERROR: GENOME %s DOES NOT HAVE FEATURE %s' % (self.name(),feature))
            sys.exit()
        else:
            return self.__featureDict[feature]

    def hasFeature(self,feature):

        if feature in self.__featureDict:
            return True
        else:
            return False

#================================================================================
#=================================FUNCTIONS======================================
#================================================================================


def loadGenome(genome_build):

    '''
    loads annotation for a genome into a genome object
    '''

    #this nested dictionary has all of the useful information and likely will have to be
    #edited so it can be configured any time
    genome_build = string.upper(genome_build)
    genomeDict = {
        'HG19':{'annot_file':'%sannotation/hg19_refseq.ucsc' % (pipeline_dir),
                'genome_directory':'/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/',
                'tf_file':'%sannotation/TFlist_NMid_hg19.txt' % (whereAmI),     
                'mask_file':'/%sgenomes/Homo_sapiens/UCSC/hg19/Annotation/Masks/hg19_encode_blacklist.bed',
                'motif_convert':'%sannotation/MotifDictionary.txt' % (whereAmI),
                #'motif_convert':'/grail/projects/medullo_final/bordo/MotifDictionary_NRL.txt',
                'motif_database':'%sannotation/VertebratePWMs.txt' % (whereAmI),
                },
        'RN6':{'annot_file':'%sannotation/rn6_refseq.ucsc' % (pipeline_dir),
                'genome_directory':'%sgenomes/Rattus_norvegicus/UCSC/rn6/Sequence/Chromosomes/' % (pipeline_dir),
                'tf_file':'%sannotation/TFlist_NMid_rn6.txt' % (whereAmI),      
                'motif_convert':'%sannotation/MotifDictionary.txt' % (whereAmI),
                'motif_database':'%sannotation/VertebratePWMs.txt' % (whereAmI),
                }

        }
    
    if genome_build not in genomeDict:
        print('ERROR: UNSUPPORTED GENOME BUILD %s. EXITING NOW' % (genome_build))
        sys.exit()

    #now attempt to load the genome
    genome = Genome(genome_build,genomeDict[genome_build]['genome_directory'],genomeDict[genome_build]['annot_file'])

    #adding additional optional features
    genome.addFeature('tf_file',genomeDict[genome_build]['tf_file'])
    if genome_build == 'HG19':
        genome.addFeature('mask',genomeDict[genome_build]['mask_file'])
    genome.addFeature('motif_convert',genomeDict[genome_build]['motif_convert'])
    genome.addFeature('motif_database',genomeDict[genome_build]['motif_database'])

    return genome



######################
# Functions
######################



def geneToEnhancerDict(genome, enhancer_file, activity_path):
    '''                                                           
    Assign each Super-Enhancer to the closest active TSS to its center
    Return a dictionary keyed by TF that points to a list of loci 
    '''
    print('Identifying enhancers and target genes from %s' %(enhancer_file))
    #should this do gene assignment????
    #for now assume gene assignment has been done
    #can later toggle to do gene assignment

    #first load the TF lists

    tf_table = utils.parseTable(genome.returnFeature('tf_file'), '\t')

    motif_table = utils.parseTable(genome.returnFeature('motif_convert'),'\t')
    #this gives all tfs that have a motif
    motif_tfs = utils.uniquify([line[1] for line in motif_table])
    #intersect w/ the activity table
    if len(activity_path) > 0:
        activity_table = utils.parseTable(activity_path,'\t')
        active_gene_list = [string.upper(line[0]) for line in activity_table]
        tf_list_refseq = [line[0] for line in tf_table if active_gene_list.count(line[1]) > 0 and motif_tfs.count(line[1]) > 0]
        tf_list_name = utils.uniquify([line[1] for line in tf_table if active_gene_list.count(line[1]) > 0 and motif_tfs.count(line[1]) > 0])
    else:
        tf_list_refseq = [line[0] for line in tf_table if motif_tfs.count(line[1]) >0]
        tf_list_name = [line[1] for line in tf_table if motif_tfs.count(line[1]) >0]
    
    print('Identified %s TFs from %s that have motifs' % (len(tf_list_name),genome.returnFeature('tf_file')))

    #keyed by gene with loci objects in the list
    gene_to_enhancer_dict = defaultdict(list)
    enhancer_to_gene_dict = defaultdict(list)

    #assuming id,chrom,start,stop w/ gene names in the last 3 columns per standard ROSE output
    enhancer_table = utils.parseTable(enhancer_file,'\t')
    print('Analyzing %s cis-regulatory regions' % (len(enhancer_table)))

    #now let's make the enhancer table by region and then by gene
    enhancerTable = [['ENHANCER_ID','CHROM','START','STOP','GENE_LIST']]
    enhancerTFTable = [['ENHANCER_ID','CHROM','START','STOP','GENE_LIST']]
    geneTable = [['GENE','TF','CHROM','START','STOP','ENHANCER_ID']]
    geneTFTable = [['GENE','CHROM','START','STOP','ENHANCER_ID']]
    geneSummaryTable = [['GENE','TF','ENHANCER_LIST']]

    #will need to track which ones are TFs
    candidate_tf_list = []
    #find the columns for gene assignment
    header = enhancer_table[0]
    header_length = len(enhancer_table[0])
    closest_index = header.index('CLOSEST_GENE')
    proximal_index = header.index('PROXIMAL_GENES')
    overlap_index = header.index('OVERLAP_GENES')

    for line in enhancer_table[1:]:
        if len(line) != header_length: #don't bother trying to figure out lines w/o target genes
            continue
        enhancer_locus = utils.Locus(line[1],line[2],line[3],'.',line[0])
        closest_gene_list = line[closest_index].split(',')
        proximal_gene_list = line[proximal_index].split(',')
        overlap_gene_list = line[overlap_index].split(',')
        all_gene_list = closest_gene_list + proximal_gene_list + overlap_gene_list
        all_gene_list = [string.upper(gene) for gene in all_gene_list]

        #gets a unique list of all tfs
        if len(activity_path) > 0:
            all_gene_list = utils.uniquify([gene for gene in all_gene_list if active_gene_list.count(gene) > 0])
        else:
            all_gene_list = utils.uniquify(all_gene_list)
        candidate_gene_list = utils.uniquify([gene for gene in all_gene_list if tf_list_name.count(gene) > 0])
        if len(all_gene_list) > 0:
            for gene in all_gene_list:

                gene_to_enhancer_dict[gene].append(enhancer_locus)
                enhancer_to_gene_dict[enhancer_locus].append(gene)
            newLine = line[0:4] + [','.join(all_gene_list)]
        else:
            newLine = line[0:4] + ['']
        enhancerTable.append(newLine)
        
        if len(candidate_gene_list) > 0:
            tfLine = line[0:4] + [','.join(candidate_gene_list)]
            enhancerTFTable.append(tfLine)

    #now iterate through each gee and list the enhancers
    gene_list = gene_to_enhancer_dict.keys()
    gene_list.sort()
    for gene in gene_list:
        if tf_list_name.count(gene) > 0:
            tf_status = 1
            candidate_tf_list.append(gene)
        else:
            tf_status = 0
        enhancer_loci = gene_to_enhancer_dict[gene]
        enhancerString =','.join([enhancer.ID() for enhancer in enhancer_loci])
        geneSummaryTable.append([gene,tf_status,enhancerString])
        for enhancer in enhancer_loci:
            newLine = [gene,tf_status,enhancer.chr(),enhancer.start(),enhancer.end(),enhancer.ID()]
            geneTable.append(newLine)
            if tf_status == 1:
                newLine = [gene,enhancer.chr(),enhancer.start(),enhancer.end(),enhancer.ID()]
                geneTFTable.append(newLine)

    return geneTable,geneTFTable,enhancerTable,enhancerTFTable,geneSummaryTable,candidate_tf_list,gene_to_enhancer_dict


def gaussianSmooth(readList, degree=5):
    '''
    Smoothing function for raw bamliquidator output
    '''

    window=degree*2-1
    weight=numpy.array([1.0]*window)
    weightGauss=[]

    for i in range(window):

        i = i-degree+1
        frac = i/float(window)
        gauss = 1/(numpy.exp((4*(frac))**2))
        weightGauss.append(gauss)

    weight=numpy.array(weightGauss)*weight
    smoothed=[0.0]*(len(readList)-window)

    for i in range(len(smoothed)):
        smoothed[i]=sum(numpy.array(readList[i:i+window])*weight)/sum(weight)

    smoothed = [0,0,0,0,0] + smoothed + [0,0,0,0] # return an array of the same length

    return smoothed


def scoreValley(locus, bam_list,max_read_length, projectName, projectFolder):
    '''
    calculate valley scores for a locus
    based on this refernce:
    http://bioinformatics.oxfordjournals.org/content/26/17/2071.full
    '''

    #takes in a bamDict where for each dataset you hold the path, the mmr, and the readlength so you can match extension
    #pull in w subprocess
    #average density for all in the group
    nBins = locus.len()/10

    #want to make an average density gram
    density_matrix = []
    for bam in bam_list:
        #calculate the extension
        extension = max_read_length - bam.getReadLengths()[0]        
        #this gives the normalized signal vector
        signal_vector = bam.liquidateLocus(locus,nBins,'.',extension,mmr = True)
        density_matrix.append(signal_vector)

    #now get the average
    if len(density_matrix) > 1:
        density = [numpy.average([line[i] for line in density_matrix]) for i in range(len(density_matrix[0]))]
    else:
        density = density_matrix[0]

    smoothDensity =  gaussianSmooth(density, 5)

    scoreArray = []
    regionMax = max(smoothDensity)

    #Now take the smooth reads and calaculate a valley score

    for i in range(len(smoothDensity)):
        score = 0
        try:
            leftmax = max(smoothDensity[i-25:i-10])
        except:
            leftmax = 'edge'
        try:
            rightmax = max(smoothDensity[i+10:i+25])
        except:
            rightmax = 'edge'

        if rightmax == 'edge' and leftmax == 'edge':
            shoulderHeightMin = 0
            shoulderHeightMax = 0
        elif leftmax == 'edge':
            shoulderHeightMin = rightmax
            shoulderHeightMax = rightmax
        elif rightmax == 'edge':
            shoulderHeightMin = leftmax
            shoulderHeightMax = leftmax
        else:
            shoulderHeightMin = min(leftmax, rightmax)
            shoulderHeightMax = max(leftmax, rightmax)

        ratio = (shoulderHeightMax-float(smoothDensity[i]))/regionMax
        if ratio > 0.3:
            score = 1
        else:
            score = 0

        scoreArray.append(score)

    return scoreArray

def stitchValleys(valleyList):
    '''
    takes a list of valley loci
    returns a stitched list of valleys to extract seq from
    '''

    valleyCollection = utils.LocusCollection(valleyList,1)
    stitchedValleyCollection = valleyCollection.stitchCollection()
    loci = []
    regions = []
    for valley in stitchedValleyCollection.getLoci():
        if [valley.chr(), valley.start(), valley.end()] not in regions:
            loci.append(valley)
            regions.append([valley.chr(), valley.start(), valley.end()])
    return loci

def findValleys(gene_to_enhancer_dict, bamFileList, projectName, projectFolder, cutoff = 0.2):
    '''
    takes in the super dict
    returns a dictionary of refseqs with all valley loci that are associated
    returns 2 kinds of bed files...
    1 = all 
    '''

    #first make the bamDict


    all_valley_bed = []
    valleyDict = {}

    #start w/ a bamFileList and make a list of bam type objects
    bam_list = [utils.Bam(bam_path) for bam_path in bamFileList]
    max_read_length = max([bam.getReadLengths()[0] for bam in bam_list])

    gene_list = gene_to_enhancer_dict.keys()
    gene_list.sort()
    ticker = 0
    print('number of regions processed:')
    for gene in gene_list:
        
        valleyDict[gene] = []

        for region in gene_to_enhancer_dict[gene]:
            if ticker %100 == 0:
                print(ticker)
            ticker+=1
            scoreArray = scoreValley(region, bam_list,max_read_length,projectName, projectFolder)
            for index,score in enumerate(scoreArray):
                if score > cutoff:
                    valley = utils.Locus(region.chr(), region.start() + index*10,
                                         region.start() + (index+1)*10, '.')
                    valleyDict[gene].append(valley)

        stitchedValleys = stitchValleys(valleyDict[gene])
        for valley in stitchedValleys:
            all_valley_bed.append([valley.chr(), valley.start(), valley.end()])
            valleyDict[gene] = stitchedValleys


    all_bed_path = projectFolder + projectName + '_all_valleys.bed'
    utils.unParseTable(all_valley_bed, all_bed_path, '\t')


    return all_bed_path





def filterSubpeaks(subpeakFile,gene_to_enhancer_dict, analysis_name,output_folder):
    '''
    takes the initial subpeaks in, stitches them, 
    '''


    # stitch the subpeaks
    print(subpeakFile)
    subpeakCollection = utils.importBoundRegion(subpeakFile,'%s_subpeak' % (analysis_name))
    
    subpeakCollection = subpeakCollection.stitchCollection()
    
    subpeakLoci = subpeakCollection.getLoci()


    all_sub_bed = []
    for locus in subpeakLoci:
        bed_line = [locus.chr(),locus.start(),locus.end(),'.',locus.ID()]
        all_sub_bed.append(bed_line)


    all_bed_path = output_folder + analysis_name + '_all_subpeak.bed'
    utils.unParseTable(all_sub_bed, all_bed_path, '\t')

    return all_bed_path



def generateSubpeakFASTA(gene_to_enhancer_dict, subpeaks, genome, projectName, projectFolder, constExtension):
    '''
    from a BED file of constituents
    generate a FASTA for the consituients contained within the canidate supers
    '''
    genomeDirectory = genome.directory()
    subpeakDict = {}
    subpeakBED = [['track name=' + projectName + ' color=204,0,204']]
    subpeakTable = utils.parseTable(subpeaks, '\t')

    subpeakLoci = [utils.Locus(l[0], int(l[1]), int(l[2]), '.') for l in subpeakTable]
    subpeakCollection = utils.LocusCollection(subpeakLoci, 50)

    for gene in gene_to_enhancer_dict.keys():
        subpeakDict[gene] = []
        for region in gene_to_enhancer_dict[gene]:
            overlaps = subpeakCollection.getOverlap(region)
            extendedOverlaps = [utils.makeSearchLocus(x, constExtension, constExtension) for x in overlaps]

            overlapCollectionTemp = utils.LocusCollection(extendedOverlaps, 50)
            overlapCollection = overlapCollectionTemp.stitchCollection()
            for overlap in overlapCollection.getLoci():
                subpeakBED.append([overlap.chr(), overlap.start(), overlap.end()])
                subpeakDict[gene].append(overlap)

    fasta = []

    for gene in subpeakDict:
        for subpeak in subpeakDict[gene]:

            fastaTitle = gene + '|'  + subpeak.chr() + '|' + str(subpeak.start()) + '|' + str(subpeak.end())
            fastaLine = utils.fetchSeq(genomeDirectory, subpeak.chr(), int(subpeak.start()+1), 
                                       int(subpeak.end()+1))

            fasta.append('>' + fastaTitle)
            fasta.append(string.upper(fastaLine))


    return subpeakBED,fasta




def makeMotifBackground(subpeakFasta,projectFolder,projectName):

    '''
    makes a 1st order markov background file for fimo
    '''

    bgCmd = 'fasta-get-markov -m 1 < ' + subpeakFasta + '  > ' + projectFolder + projectName + '_bg.meme'
    bg_path = '%s%s_bg.meme' %(projectFolder,projectName)
    subprocess.call(bgCmd, shell=True)

    return bg_path


def findMotifs(subpeakFasta,bg_path,candidate_tf_list, projectFolder, analysis_name, motifConvertFile, motifDatabaseFile):
    '''
    takes the refseq to subpeak seq dict
    returns the networkx object with all connections
    '''
    fimoFolder = utils.formatFolder(projectFolder + 'FIMO/', True)
    subpeak_name = subpeakFasta.split('/')[-1].split('.')[0]
    output = '%s%s_fimo.txt'  % (fimoFolder,subpeak_name)
    # Create a dictionary to call motif names keyed on gene names
    motifDatabase = utils.parseTable(motifConvertFile, '\t')
    motifDatabaseDict = {} #create a dict keyed by TF with multiple motifs

    for line in motifDatabase:
        motifDatabaseDict[line[1]] = []
    for line in motifDatabase:
        motifDatabaseDict[line[1]].append(line[0])

    candidate_tf_list.sort()

    #now make a list of all motifs
    motif_list = []
    for tf in candidate_tf_list:
        motif_list += motifDatabaseDict[tf]

    motif_list = utils.uniquify(motif_list)

    fimo_bash_path = '%s%s_fimo.sh' % (fimoFolder,analysis_name)
    fimo_bash = open(fimo_bash_path,'w')
    fimo_bash.write('#!/usr/bin/bash\n\n')

    fimoCmd = 'fimo'
    for motif in motif_list:
        fimoCmd += ' --motif ' + "'%s'" % (str(motif))


    #fimoCmd += ' --thresh 1e-5' #if you want to increase stringency
    fimoCmd += ' -verbosity 1'  # thanks for that ;)!
    fimoCmd += ' -text'
    fimoCmd += ' -oc ' + projectFolder + 'FIMO'
    fimoCmd += ' --bgfile %s' % (bg_path)
    fimoCmd += ' ' + motifDatabaseFile + ' '
    fimoCmd += subpeakFasta
    fimoCmd += ' > '+ output
    print fimoCmd
    fimo_bash.write(fimoCmd)
    fimo_bash.close()

    fimoOutput = subprocess.call(fimoCmd, shell=True)  #will wait that fimo is done to go on

    return output


#after fimo, we should collapse the motif edges into beds before going into networks

def collapseFimo(fimo_output,gene_to_enhancer_dict,candidate_tf_list,output_folder,analysis_name,motifConvertFile):

    '''
    collapses motifs from fimo
    for each source node (TF) and each target node (gene enhancer regions), collapse motif instances
    then spit out a ginormous set of beds and a single crazy collapsed bed
    '''
    
    #first build up the motif name conversion database

    motifDatabase = utils.parseTable(motifConvertFile, '\t')
    motifDatabaseDict = defaultdict(list)
    # The reverse of the other dict, from motif name to gene name
    # a motif can go to multiple genes
    for line in motifDatabase:
        motifDatabaseDict[line[0]].append(line[1])



    #make the folder to store motif beds
    utils.formatFolder('%smotif_beds/' % (output_folder),True)

    edgeDict = {}
    #first layer are source nodes
    for tf in candidate_tf_list:
        edgeDict[tf] = defaultdict(list) #next layer are target nodes which are derived from the fimo output
        

    fimoTable = utils.parseTable(fimo_output,'\t')

    for line in fimoTable[1:]:
        source_tfs = motifDatabaseDict[line[0]]   #motifId
        for source in source_tfs:
            if candidate_tf_list.count(source) == 0:
                continue
            region = line[1].split('|')

            target = region[0]
            target_locus = utils.Locus(region[1],int(region[2]) + int(line[2]), int(region[2]) + int(line[3]),'.')
        
            #what's missing here is the enhancer id of the target locus
            try:
                edgeDict[source][target].append(target_locus)
            except KeyError:
                print('this motif is not in the network')
                print(line)
                sys.exit()


    #now we actually want to collapse this down in a meaningful way
    #overlapping motifs count as a single binding site. This way a TF with tons of motifs
    #that finds the same site over and over again doesn't get over counted
    all_bed = []
    all_bed_path = '%s%s_all_motifs.bed' % (output_folder,analysis_name)
    for tf in candidate_tf_list:
        print(tf)
        target_nodes = edgeDict[tf].keys()
        bed_header = ['track name = "%s" description="%s motifs in %s"' % (tf,tf,analysis_name)]
        all_bed.append(bed_header)
        target_bed = [bed_header]
        target_bed_path = '%smotif_beds/%s_motifs.bed' % (output_folder,tf)
        for target in target_nodes:
            edgeCollection = utils.LocusCollection(edgeDict[tf][target],50)
            edgeCollection = edgeCollection.stitchCollection()
            edgeLoci = edgeCollection.getLoci()
            edgeDict[tf][target] = edgeLoci
            for locus in edgeLoci:
                bed_line = [locus.chr(),locus.start(),locus.end(),target,'','+']
                target_bed.append(bed_line)
                all_bed.append(bed_line)

        utils.unParseTable(target_bed,target_bed_path,'\t')

    #now the loci are all stitched up 
    utils.unParseTable(all_bed,all_bed_path,'\t')
    return edgeDict
                                   


def buildGraph(edgeDict,gene_to_enhancer_dict,output_folder, analysis_name,cutoff=1):
    '''
    from the collapsed edge dictionary, build a target graph
    require at least n motifs to constitute an edge where n is set by cutoff. 
    default is 1
    '''

    node_list = edgeDict.keys()
    node_list.sort()
    #this is only edges between TFs
    graph = nx.DiGraph(name=analysis_name)
    graph.add_nodes_from(node_list)
    

    #this stores ALL edges identified by motifs
    edge_table = [['SOURCE','TARGET','CHROM','START','STOP','REGION_ID','TF_INTERACTION']]
    edge_output = '%s%s_EDGE_TABLE.txt' % (output_folder,analysis_name)

    for source in node_list:
        print(source)
        target_list = edgeDict[source].keys()
        target_list.sort()
        for target in target_list:

            #now we need to see which target regions this guy overlaps
            target_regions = gene_to_enhancer_dict[target]
            target_collection = utils.LocusCollection(target_regions,50)

            #get the edges hitting that target
            edgeLoci = edgeDict[source][target]
            if node_list.count(target) > 0:
                tf_interaction = 1
            else:
                tf_interaction = 0
            #only add to the graph if this is a TF/TF interaction
            if len(edgeLoci) >= cutoff and node_list.count(target) > 0:
                graph.add_edge(source,target)
                
            #now for each edge, add to the table
            for edgeLocus in edgeLoci:
                regionString = ','.join([locus.ID() for locus in target_collection.getOverlap(edgeLocus)])
                edgeLine = [source,target,edgeLocus.chr(),edgeLocus.start(),edgeLocus.end(),regionString,tf_interaction]
                edge_table.append(edgeLine)

    utils.unParseTable(edge_table,edge_output,'\t')
    return graph






def formatNetworkOutput(graph, output_folder, analysis_name, candidate_tf_list):
    '''
    takes the networkx graph
    returns all figures, tables, etc
    '''

    # output the network as a .ntx dictionary of lists

    networkFilename = output_folder + analysis_name + '.ntx'
    networkFile = open(networkFilename, 'w')
    networkDictOfLists = nx.to_dict_of_lists(graph)
    pickle.dump(networkDictOfLists, networkFile)

    # output the adjacency list and nodelist
    nodeFile = output_folder + analysis_name + '_NODELIST.txt'
    nodeList = [ [n] for n in graph.nodes_iter()]
    utils.unParseTable(nodeList, nodeFile, '\t')

    adjFile = output_folder + analysis_name + '_ADJ_LIST.txt'
    adjList = graph.adjacency_list()
    utils.unParseTable(adjList, adjFile, '\t')

    edgesTable = [['From', 'To']]
    targetList = []
    for i,gene in enumerate(nodeList):
        for j in adjList[i]:
            newline = [gene[0],j]
            edgesTable.append(newline)
            TFname = gene[0]

    edgeFile = output_folder + analysis_name + '_EDGE_LIST.txt'
    utils.unParseTable(edgesTable, edgeFile, '\t')


    # Make the degree table    
    degTable = [['Tf', 'In_Degree', 'Out_Degree', 'Total_Connections' ]]
    degFile = output_folder + analysis_name + '_DEGREE_TABLE.txt'

    for node in graph.nodes(): #shouldn't we output the table for the TFs that have motifs only ? for canidateMotifs in graph.nodes()....
        newline = [node, graph.in_degree()[node], graph.out_degree()[node], graph.degree()[node]]
        degTable.append(newline)

    utils.unParseTable(degTable, degFile, '\t')

    print 'DEFINING THE CORE REGULATORY CIRCUIT'

    autoreg = graph.selfloop_edges()
    selfLoops = [x for x,y in autoreg]
    selfLoopFile = output_folder + analysis_name + '_SELF_LOOPS.txt'
    utils.unParseTable(selfLoops, selfLoopFile, '')

    #recover bidirectional edges

    pairs = []
    for n in selfLoops:
        for m in selfLoops:
            if n != m:
                if graph.has_edge(n,m) and graph.has_edge(m,n):
                    pairs.append([n,m])
    
    unDirGraph = nx.from_edgelist(pairs)
    cliqueGen = find_cliques_recursive(unDirGraph)
    cliqueList = list(cliqueGen)

    utils.unParseTable(cliqueList, output_folder + analysis_name + '_CLIQUES_ALL.txt', '\t')

    cliqueRanking = []
    outDegreeDict = graph.out_degree()

    for c in cliqueList:
        score = 0
        for gene in c:
            score += outDegreeDict[gene]
        score = score/len(c)
        if score > 0 and len(c) > 2:
            cliqueRanking.append((c, score))


    sortCliqueRanking = sorted(cliqueRanking, reverse=True, key=lambda x:x[1])
    cliqueFile = output_folder + analysis_name + '_CLIQUE_SCORES_DEGREE.txt'
    utils.unParseTable(sortCliqueRanking, cliqueFile, '\t')

    factorEnrichmentDict = {}

    for factor in selfLoops:
        factorEnrichmentDict[factor] = 0
    for pair in cliqueRanking:
        c = pair[0]
        for factor in c:
            factorEnrichmentDict[factor] += 1

    factorRankingTable = []
    for factor in selfLoops:
        newline = [factor, factorEnrichmentDict[factor]/float(len(cliqueRanking))]
        factorRankingTable.append(newline)

    factorRankingFile = output_folder + analysis_name + '_ENRICHED_CLIQUE_FACTORS.txt'
    utils.unParseTable(factorRankingTable, factorRankingFile, '\t')

    # Begin VSA scoring 

    # Initiate the graph
    G=nx.Graph()

    #recover bidirectional edges
    bidirectionalEdges = pairs

    #fill up the graph
    G.add_nodes_from(selfLoops)
    G.add_edges_from(bidirectionalEdges)

    #find all the cliques
    cliques = find_cliques_recursive(G)
    cliqueList = list(cliques)

    print 'Number of cliques:'
    print len(cliqueList)

    #count the occurences of the TFs accross the loops

    dicoTFinloopsCounts={}

    for clique in cliqueList:
        for TF in clique:

            if dicoTFinloopsCounts.has_key(TF):
                dicoTFinloopsCounts[TF]+=1

            else:
                dicoTFinloopsCounts[TF]=1

    #calculate a score by loop

    cliqueRanking = []

    cliqueNub = 0


    for clique in cliqueList:
        cliqueScore=0


        for TF in clique:
            cliqueScore = (float(cliqueScore) + (float(dicoTFinloopsCounts[TF])))
            cliqueRanking.append((clique, cliqueScore/len(clique), len(clique)))

    #print(cliqueRanking)
    sortCliqueRanking = sorted(cliqueRanking, reverse=True, key=lambda x:x[1])
    #print(sortCliqueRanking)
    cliqueFile = output_folder + analysis_name + '_CLIQUE_SCORES_VSA.txt'
    utils.unParseTable(sortCliqueRanking, cliqueFile, '\t')

    print 'Top CRC:'
    print sortCliqueRanking[0]

    # Visualizations

    #sizeFile = output_folder + analysis_name + '_CANIDATE_TF_AND_SUPER_TABLE.txt'
    #os.system('Rscript networkScatter.R ' + degFile + ' ' + sizeFile + ' ' +
    #          output_folder + analysis_name + '_NETWORK_SCATTER.pdf')



###################### 
#
# Main Method
#
######################

def main():

    import argparse
    parser = argparse.ArgumentParser(usage="usage: prog [options] -e [ENHANCER_FILE] -b [BAM_FILE] -g [GENOME] -o [OUTPUTFOLDER] -n [NAME]" )



    #required flags                                                                                                                   
    parser.add_argument("-e","--enhancer_file", dest="enhancers", default=None,type=str,
                        help = "Provide a ROSE generated enhancer table (_AllEnhancers.table.txt)",required=True)

    parser.add_argument("-g","--genome",dest="genome", default = None,type=str,
                        help = "Provide the build of the genome to be used for the analysis. Currently supports HG19, HG18 and MM9",required=True)
    parser.add_argument("-o","--output",dest="output", default = None,type=str,
                        help = "Enter an output folder",required=True)
    parser.add_argument("-n","--name",dest="name", default = None,type=str,
                        help = "Provide a name for the job",required=True)


    #you either need bams for valleys or subpeaks
    parser.add_argument("-b","--bam",dest="bam", default = None,type=str,
                        help = "Enter a comma separated list of bams of valley finding",required=False)
    parser.add_argument("-s","--subpeaks", dest="subpeaks",default=None,type=str,
                        help = "Enter a BED file of regions to search for motifs",required=False)



    #additional options                                                                                  
    parser.add_argument("-a","--activity",dest="activity", default = None,type=str,
                        help = "A table with active gene names in the first column",required=False)
    parser.add_argument("-l","--extension-length", dest="extension", default=100,type=int,
                        help = "Enter the length to extend subpeak regions for motif finding",required=False)
    parser.add_argument("-B","--background", dest="background", default=None,type=str,
                        help = "Provide a background BAM file",required=False)
    parser.add_argument("-N", "--number", dest="number", default=1,type=int,
                        help = "Enter the number of non overlapping motifs in a region required to assign a binding event. Default=1",required=False)     #I have modified the destination of -N option so that it is different from the destination of -E option
    parser.add_argument("--motifs", dest="motifs", default=False,type=str,
                        help = "Enter additional PWM file for the analysis",required=False)
    parser.add_argument("-t","--tfs", dest="tfs",default=None,type=str,
                        help = "Enter additional TFs (comma separated) to be used in the bindinf analysis",required=False)

    args = parser.parse_args()





    #=====================================================================================
    #===============================I. PARSING ARGUMENTS==================================
    #=====================================================================================


    ###
    # Define all global file names
    ###
    print(args)
    genome = loadGenome(args.genome)

    motifDatabaseFile = genome.returnFeature('motif_database')
    motifConvertFile = genome.returnFeature('motif_convert')

    # User input files
    enhancer_file = args.enhancers

    if args.bam == None and args.subpeaks == None:
        print('ERROR: Must provide either bams for valley finding or subpeaks as a .bed')
        sys.exit()

    #set the subpeak file
    if args.subpeaks:
        subpeakFile = args.subpeaks
    else: subpeakFile = None


    #will need to fix bams down the line to take in multiple bams
    if args.bam:
        bamFileList = [bam_path for bam_path in args.bam.split(',') if len(bam_path) >0]
        print(bamFileList)
    else:
        bamFileList = []

    if args.background:
        background = args.background

    else: 
        background = None


    #output folder and analysis name
    output_folder = utils.formatFolder(args.output,True)
    analysis_name = args.name

    #optional arguments
    #activity path
    activity_path = args.activity

    #motif extension
    constExtension = args.extension

    print('\n\n#======================================\n#===========I. DATA SUMMARY============\n#======================================\n')

    print('Analyzing TF connectivity for %s' % (analysis_name))
    print('Writing output to %s' % (output_folder))
    if subpeakFile:
        print('Using %s to define subpeaks for motif finding' % (subpeakFile))
    else:
        print('Identifying valleys from .bam files')
    print('Using %s to define active genes' % (activity_path))


    #=====================================================================================
    #=======================II. IDENTIFYING CANDIDATE TFS AND NODES=======================
    #=====================================================================================

    print('\n\n#======================================\n#===II. MAPPING GENES AND ENHANCERS====\n#======================================\n')
    
    geneTable,geneTFTable,enhancerTable,enhancerTFTable,geneSummaryTable,candidate_tf_list,gene_to_enhancer_dict= geneToEnhancerDict(genome, enhancer_file, activity_path)
    #write these guys to disk

    gene_out = '%s%s_GENE_TABLE.txt' % (output_folder,analysis_name)
    gene_tf_out = '%s%s_GENE_TF_TABLE.txt' % (output_folder,analysis_name)

    enhancer_out = '%s%s_ENHANCER_TABLE.txt' % (output_folder,analysis_name)
    enhancer_tf_out = '%s%s_ENHANCER_TF_TABLE.txt' % (output_folder,analysis_name)

    summary_out= '%s%s_GENE_SUMMARY.txt' % (output_folder,analysis_name)
    
    utils.unParseTable(enhancerTable,enhancer_out,'\t')    
    utils.unParseTable(enhancerTFTable,enhancer_tf_out,'\t')

    utils.unParseTable(geneTable,gene_out,'\t')
    utils.unParseTable(geneTFTable,gene_tf_out,'\t')

    utils.unParseTable(geneSummaryTable,summary_out,'\t')
    

    print('Identified %s genes w/ proximal cis-regulatory elements' % (len(gene_to_enhancer_dict)))
            
    print('Identified %s candidate TFs' % (len(candidate_tf_list)))
    print(candidate_tf_list)


    #=====================================================================================
    #==========================III. FINDING VALLEYS/SUBPEAKS==============================
    #=====================================================================================

    print('\n\n#======================================\n#=====III. FINDING VALLEYS/SUBPEAKS====\n#======================================\n')


    #so here we would need to find valleys everywhere
    if subpeakFile == None:
        print('finding valleys')
        #note: the tf_bed_path is for networks, all is for out degree finding
        all_bed_path = findValleys(gene_to_enhancer_dict, bamFileList, analysis_name, output_folder, cutoff = 0.2)
    else:
        print('Using subpeaks from %s' % (subpeakFile))
        all_bed_path = filterSubpeaks(subpeakFile,gene_to_enhancer_dict,analysis_name,output_folder)


    #first make the subpeak bed and subpeak fasta for the tfs
    all_sub_bed,all_fasta = generateSubpeakFASTA(gene_to_enhancer_dict, all_bed_path, genome, analysis_name,output_folder, constExtension)

    if subpeakFile == None:
        #this is the case where we did valleys #only reason you would need to output the sub bed
        all_sub_out = '%s%s_all_subpeak.bed' % (output_folder,analysis_name)
        utils.unParseTable(all_sub_bed,all_sub_out,'\t')


    #writing the all subpeak fasta out to disk
    all_fasta_out = '%s%s_all_subpeak.fasta' % (output_folder,analysis_name)
    utils.unParseTable(all_fasta,all_fasta_out,'')
        

    #=====================================================================================
    #=================================IV. FINDING MOTIFS==================================
    #=====================================================================================

    print('\n\n#======================================\n#======IV. RUNNING MOTIF FINDING=======\n#======================================\n')


    #first make background
    bg_path = makeMotifBackground(all_fasta_out,output_folder,analysis_name)

    #find motifs for all regions
    fimo_out = findMotifs(all_fasta_out,bg_path,candidate_tf_list, output_folder, analysis_name, motifConvertFile, motifDatabaseFile)

    edgeDict = collapseFimo(fimo_out,gene_to_enhancer_dict,candidate_tf_list,output_folder,analysis_name,motifConvertFile)


    #=====================================================================================
    #============================V. RUNNING NETWORK ANALYSIS==============================
    #=====================================================================================

    print('\n\n#======================================\n#========V. BUILDING NETWORK============\n#======================================\n')


    print('building graph and edge table')
    graph = buildGraph(edgeDict,gene_to_enhancer_dict,output_folder, analysis_name,cutoff=1)

    formatNetworkOutput(graph, output_folder, analysis_name, candidate_tf_list)

        
    print('yay')

    sys.exit()




if __name__ == '__main__':
    main()

