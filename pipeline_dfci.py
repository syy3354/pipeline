#!/usr/bin/env python
#pipeline.py


'''
The MIT License (MIT)

Copyright (c) 2013 Charles Lin

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

#module of functions/code/structures from the myc project that has now been
#addapted for general use



#==========================================================================
#==========================DEPENDENCIES====================================
#==========================================================================




import sys
import os
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipelineFolder = '%s/' % (whereAmI) # need to set this to where this code is stored

sys.path.append(pipelineFolder)

print('\nUsing following version of python:\n')
print(sys.version)
print('\n\n')



from utils import *
import utils
import datetime
import subprocess
import time
import re
import random
import string
import numpy


from collections import defaultdict

#==========================================================================
#==============================GLOBAL PATHS================================
#==========================================================================


samtoolsString = 'samtools'
bamliquidator_path = 'bamliquidator_batch.py'

fastqDelimiter = '::' #delimiter for pairs in fastqs




#==========================================================================
#===========================TABLE OF CONTENTS==============================
#==========================================================================

#CODE IN THIS MODULE IS SPLIT UP INTO SEVERAL SECTIONS


#-------------------------------------------------------------------------#
#                                                                         #
#                         FORMATTING AND DATA INPUT                       #
#                                                                         #
#-------------------------------------------------------------------------#


#FORMATTING FOLDERS
#formatFolder(folderName,create=False)

#FORMATTING THE MASTER DATA TABLE
#formatDataTable(dataFile):

#FORMATTING FUNCTIONS
#makePipelineTable(sampleTableFile,dirPath,bamPath,outputFile,overwrite=False):
#makeGenialisTable(dataFile,outFilePath,organism='',seqType='',paired=False,collection='',annotator='',source='',strain='',tissue='',age='',genotype='',molecule='',libraryStrategy='',exPro='',libraryConst='',other1='',other2='')

#LOADING THE MASTER DATA TABLE
#def loadDataTable(dataFile):
#def writeDataTable(dataDict,outFile):
#def summary(dataFile,outputFile=''):
#def makeBamTable(dataFile,output):


#-------------------------------------------------------------------------#
#                                                                         #
#                          TONY DATABASE TOOLS                            #
#                                                                         #
#-------------------------------------------------------------------------#

#INTERACTING WITH TONY
#getTONYInfo(uniqueID,column =''):


#-------------------------------------------------------------------------#
#                                                                         #
#                                ALIGNMENT                                #
#                                                                         #
#-------------------------------------------------------------------------#


#CALLING BOWTIE TO MAP DATA
#def makeBowtieBashJobs(pipelineFile,namesList = [],launch=True,overwrite=False):
#def makeBowtieBashJobsSlurm(pipelineFile,namesList = [],launch=True,overwrite=False):
#def callBowtie(dataFile,dataList = [],overwrite = False):

#GETTING MAPPING STATS
#def bowtieStats(dataFile,namesList=[]):

#MERGING BAMS
#def mergeBams(dataFile,mergedName,namesList,color='',background =''):

#FILTERING BAMS
#def filterBams(dataFile,namesList = [],tempFolder = '/raider/BOWTIE_TEMP/',bamFolder='/grail/bam/filtered/'):

#-------------------------------------------------------------------------#
#                                                                         #
#                              PEAK FINDING                               #
#                                                                         #
#-------------------------------------------------------------------------#


#CALLING MACS
#def callMacsQsub(dataFile,macsFolder,namesList = [],overwrite=False,pvalue='1e-9'):
#def callMacs(dataFile,macsFolder,namesList = [],overwrite=False,pvalue='1e-9',useBackground =True):
#def callMacsSlurm(dataFile,macsFolder,namesList = [],overwrite=False,pvalue='1e-9',useBackground =True):
#def callMacs2(dataFile,macsFolder,namesList = [],broad=True,noBackground = False,pairedEnd = False,overwrite=False,pvalue='1e-9'):

#FORMATTING MACS OUTPUT
#def formatMacsOutput(dataFile,macsFolder,macsEnrichedFolder,wiggleFolder,wigLink='',useBackground=True):


#-------------------------------------------------------------------------#
#                                                                         #
#                              GFF TOOLS                                  #
#                                                                         #
#-------------------------------------------------------------------------#


#MAKING GFFS OF TSS REGIONS
#def makeGeneGFFs(annotFile,gffFolder,species='HG18'):

#MAKING GFFS OF CHROMS
#def makeChromGFFs(chromLengthFile,gffFolder,chromList = [],genome='HG18',binSize = 100000,singleGFF = True):

#MAKING GFFS OF ENHANCER REGIONS
#def makeEnhancerGFFs(dataFile,gffName,namesList,annotFile,gffFolder,enrichedFolder,window=2000,macs=True):

#MAKING GFFS OF ENRICHED REGIONS
#def makeEnrichedGFFs(dataFile,namesList,gffFolder,enrichedFolder,macs=True,window=0):

#MAKING GFFS OF PROMOTER REGIONS
#def makePromoterGFF(dataFile,annotFile,promoterFactor,enrichedFolder,gffFolder,window=0,transcribedGeneFile=''):


#-------------------------------------------------------------------------#
#                                                                         #
#                             MAPPING TOOLS                               #
#                                                                         #
#-------------------------------------------------------------------------#


#MAP ENRICHED REGIONS TO GFF
#def mapEnrichedToGFF(dataFile,setName,gffList,cellTypeList,enrichedFolder,mappedFolder,macs=True,namesList=[],useBackground=True):


#MAPPING BAMS TO GFFS
#def mapBams(dataFile,cellTypeList,gffList,mappedFolder,nBin = 200,overWrite =False,nameList = []):
#def mapBamsQsub(dataFile,cellTypeList,gffList,mappedFolder,nBin = 200,overWrite =False,nameList = []):
#def mapBamsBatch(dataFile,gffList,mappedFolder,overWrite =False,namesList = [],extension=200)

#FORMATTING MAPPING SIGNAL
#def makeSignalTable(dataFile,gffFile,mappedFolder,namesList = [],medianNorm=False,output =''):

#WRAPPING MAPPER
#def map_regions(dataFile,gffList,mappedFolder,signalFolder,names_list=[],medianNorm=False,output=''):
#MAKING GFF LISTS
#def makeGFFListFile(mappedEnrichedFile,setList,output,annotFile=''):

#-------------------------------------------------------------------------#
#                                                                         #
#                            PLOTTING TOOLS                               #
#                                                                         #
#-------------------------------------------------------------------------#


#PLOTTING INDIVIDUAL GENES
#def callGenePlot(dataFile,geneID,plotName,annotFile,namesList,outputFolder,region='TXN',yScale = 'UNIFORM'):

#BATCH PLOTTING REGIONS
#def callBatchPlot(dataFile,inputFile,plotName,outputFolder,namesList=[],uniform=True,bed ='',plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '',scaleFactorString =''):

#PLOTTING TF CORR HEATMAPS FROM CRC MOTIF BEDS
#def plotCRCCorrMaps(analysis_name,motifBedDir,tf_list_path='',window=50)

#-------------------------------------------------------------------------#
#                                                                         #
#                              META TOOLS                                 #
#                                                                         #
#-------------------------------------------------------------------------#



#MAKING META GFFS OF TXN REGIONS
#def makeMetaGFFs(annotFile,gffFolder,genome,geneListFile =''):

#MAPPING BAMS FOR METAS
#def mapMetaBams(dataFile,metaName,gffList,cellTypeList,metaFolder,nameList= [],overwrite=False):

#FINISHING METAS
#def finishMetas(metaFolder,settingsFileList=[]):

#MAKING ORDERED HEATMAPS
#def callHeatPlotOrdered(dataFile,gffFile,namesList,orderByName,geneListFile,outputFolder,mappedFolder,relative=False):


#-------------------------------------------------------------------------#
#                                                                         #
#                              ROSE TOOLS                                 #
#                                                                         #
#-------------------------------------------------------------------------#




#CALLING ROSE
#def callRose(dataFile,macsEnrichedFolder,parentFolder,namesList=[],extraMap = [],inputFile='',tss=2500,stitch=12500,bashFileName ='',mask=''):
#def callRose2(dataFile,macsEnrichedFolder,parentFolder,namesList=[],extraMap = [],inputFile='',tss=2500,stitch='',bashFileName ='',mask=''):
#def callRose2Slurm(dataFile,macsEnrichedFolder,parentFolder,namesList=[],extraMap = [],inputFile='',tss=2500,stitch='',bashFileName ='',mask=''):


#-------------------------------------------------------------------------#
#                                                                         #
#                               CRC TOOLS                                 #
#                                                                         #
#-------------------------------------------------------------------------#


#def call_crc(analysis_name,enhancer_path,subpeak_path,activity_path,crc_folder=''):

#-------------------------------------------------------------------------#
#                                                                         #
#                          EXPRESSION TOOLS                               #
#                                                                         #
#-------------------------------------------------------------------------#

#MAKING EXPRESSION TABLES
#def mapHisat(dataFile,namesList=[],useSRA=False,pCount=16,Launch=True):
#def makeCuffTable(dataFile,analysisName,gtfFile,cufflinksFolder,groupList=[],bashFileName = ''):
#def makeCuffTableSlurm(dataFile,analysisName,gtfFile,cufflinksFolder,groupList=[],bashFileName = ''):


#-------------------------------------------------------------------------#
#                                                                         #
#                              GECKO TOOLS                                #
#                                                                         #
#-------------------------------------------------------------------------#


#-------------------------------------------------------------------------#
#                                                                         #
#                               UCSC TOOLS                                #
#                                                                         #
#-------------------------------------------------------------------------#

#def makeTrackHub(analysis_name,project_folder,chrom_sizes, dataFileList=[], wiggle_dir='',web_dir='/storage/cylin/web/Lin_Lab_Track_Hubs/',hub_name='',hub_short_lab='',hub_long_lab='',EMAIL='',fileType='bigWig',col='0,0,0',scaled=False)



#-------------------------------------------------------------------------#
#                                                                         #
#                              GSEA TOOLS                                 #
#                                                                         #
#-------------------------------------------------------------------------#

#def wrapGSEA(gctPath,clsPath,sample_1,sample_2,analysis_name,output_folder,metric,gmxPath='',gseaPath='',launch=True)

#============================================================================================================
#============================================================================================================
#============================================================================================================




#-------------------------------------------------------------------------#
#                                                                         #
#                         FORMATTING AND DATA INPUT                       #
#                                                                         #
#-------------------------------------------------------------------------#


#==========================================================================
#==========================DATA TABLE FORMAT===============================
#==========================================================================

#the master data table
#everything starts with this
#format is as follows

#FILE_PATH UNIQUE_ID NAME BACKGROUND ENRICHED_REGION ENRICHED_MACS COLOR

#FILE_PATH = FOLDER WHERE THE BAM FILE LIVES
#UNIQUE_ID = TONY UNIQUE ID
#GENOME = GENOME USED. MM9, HG18 ARE SUPPORTED
#NAME = IDENTIFIER OF THE DATASET NEEDS TO BE UNIQUE
#BACKGROUND = NAME OF THE BACKGROUND DATASET
#ENRICHED_REGION = NAME OF THE ERROR MODEL ENRICHED REGION OUTPUT
#ENRICHED_MACS = NAME OF THE MACS PEAKS BED FILE
#COLOR = COMMA SEPARATED RGB CODE

#==========================================================================
#===================FORMATTING FOLDERS=====================================
#==========================================================================

def formatFolder(folderName,create=False):

    '''
    makes sure a folder exists and if not makes it
    returns a bool for folder
    '''
    
    if folderName[-1] != '/':
        folderName +='/'

    try: 
        foo = os.listdir(folderName)
        return folderName
    except OSError:
        print('folder %s does not exist' % (folderName))
        if create:
            os.system('mkdir %s' % (folderName))
            return folderName
        else:
                    
            return False 


#==========================================================================
#===================FORMATTING THE MASTER DATA TABLE=======================
#==========================================================================

def formatDataTable(dataFile):
    '''
    formats the dataFile and rewrites.
    first 3 columns are required for every line.
    if they aren't there the line is deleted
    '''
    print('reformatting data table')


    dataTable = parseTable(dataFile,'\t')

    newDataTable = [['FILE_PATH','UNIQUE_ID','GENOME','NAME','BACKGROUND','ENRICHED_REGION','ENRICHED_MACS','COLOR','FASTQ_FILE']]
    #first check to make sure the table is formatted correctly
    
    for line in dataTable[1:]:
        if len(line) < 3:
            continue
        #this spots header lines that may be out of place
        if line[0] == 'FILE_PATH':
            continue 
        #check if it at least has the first 3 columns filled in
        if len(line[0]) == 0 or len(line[1]) == 0 or len(line[2]) == 0:
            print('ERROR required fields missing in line')
            print(line)
        #if the first three are filled in, check to make sure there are 8 columns
        else:
            if len(line) > 3 and len(line) < 9:
                newLine = line + (8-len(line))*[''] + ['NA']
                newDataTable.append(newLine)
            elif len(line) >= 9:
                newLine = line[0:9]
                newDataTable.append(newLine)
    
    #lower case all of the genomes
    #make the color 0,0,0 for blank lines and strip out any " marks
    for i in range(1,len(newDataTable)):
        newDataTable[i][2] = string.lower(newDataTable[i][2])
        color = newDataTable[i][7]
        if len(color) == 0:
            newDataTable[i][7] = '0,0,0'
    unParseTable(newDataTable,dataFile,'\t')
    return newDataTable


#========================================================================
#=======================FORAMTTING FUNCTIONS=============================
#========================================================================


def makePipelineTable(sampleTableFile,dirPath,bamPath,outputFile,overwrite=False):

    '''
    makes a standard pipeline table in the same directory as the sample Table File
    which should be the project directory file
    uses a standard WI annotation xls
    '''
    if sampleTableFile.split('.')[-1] == 'xls':
        sampleTable = parseTable(sampleTableFile,'\t',excel=True)
    else:
        sampleTable = parseTable(sampleTableFile,'\t',excel=False)

    #check if the outputfile exists
    #if it does, append
    if overwrite:
        pipelineTable = [['FILE_PATH','UNIQUE_ID','GENOME','NAME','BACKGROUND','ENRICHED_REGIONS','ENRICHED_MACS','COLOR','FASTQ_FILE']]    
    else:
        try:
            pipelineTable = parseTable(outputFile,'\t')
        
        except IOError:
            pipelineTable = [['FILE_PATH','UNIQUE_ID','GENOME','NAME','BACKGROUND','ENRICHED_REGIONS','ENRICHED_MACS','COLOR','FASTQ_FILE']]    

    for line in sampleTable[1:]:

        if line[0] == '':
            break
        uniqueID = "%s_%s" % (line[8],line[4])
        genome = line[14]
        name = line[4]
        cellLine = string.lower(line[4].split('_')[0])
        filePath = "%s%s/%s/" % (bamPath,genome,cellLine)
        barcode = line[5].split('-')[1]
        fastqFile = "%s%s-s_%s_1_sequence.txt.tar.gz" % (dirPath,barcode,line[0])

        if name.count('WCE') == 0:
            newLine = [filePath,uniqueID,genome,name,'NONE','NONE','','0,0,0',fastqFile]

        background = line[3]
        if len(background) == 0:
            newLine = [filePath,uniqueID,genome,name,'','','','0,0,0',fastqFile]
        else:
            newLine = [filePath,uniqueID,genome,name,line[3],'','','0,0,0',fastqFile]
            
        pipelineTable.append(newLine)

    unParseTable(pipelineTable,outputFile,'\t')

def makeGenialisTable(dataFile,outFilePath,organism='',seqType='',paired=False,collection='',annotator='',source='',strain='',tissue='',age='',genotype='',molecule='',libraryStrategy='',exPro='',libraryConst='',other1='',other2=''):
    #Most of these inputs are left as blanks, because they are dependent on how the data was prepared.
    #   dataFile is the pipeline formated table used for the Lin Lab pipeline
    #   outFilePath designates where you would like the annotation table to be written to
    #   organism type is based on the genome used, but must be written with the proper format
    #   seqType is the type of sequencing (e.g. RNAseq, CHiPseq, etc)
    #   paired is default False, if you have paired end data, set this to True
    #   collection is the name of your project, it is helpful to you to make it meaningful
    #   annotator is the person preparing the annotation    
    #   source, strain, tissue, age, and genotype are self-explanitory
    #   molecule can be one of the following:
    #                  total RNA
    #                  polyA RNA
    #                  cytoplasmic RNA
    #                  nuclear RNA
    #                  genomic DNA
    #                  protein
    #                  other
    #  libraryStrategy can be used to denote strand i.e. single end ChIP-Seq
    #  exPro is the extraction protocol used i.e. young lab chip protocol
    #  libraryConst is the library construction protocol used i.e. rubicon thruplex
    #  other1 and other2 can be used to denote any other information that is not in the table
    

    #This is the header for genialis annotation tables
    genialis_header = ['NAME','FASTQ_R1','FASTQ_R2','SEQ_TYPE','PAIRED','COLLECTION','ANNOTATOR','SOURCE','ORGANISM','STRAIN','TISSUE','AGE','GENOTYPE','MOLECULE',
                       'LIBRARY_STRATEGY','EXTRACTION_PROTOCOL','LIBRARY_CONSTRUCTION_PROTOCOL','OTHER_CHAR_1','OTHER_CHAR_2']
    
    #This selects the proper organism based on the genome you are using
    genomeDict = {'HG38':'Homo sapiens',
                  'HG19':'Homo sapiens',
                  'HG19_ERCC':'Homo sapiens',
                  'MM9':'Mus musculus',
                  'MM10':'Mus musculus',
                  'RN6':'Rattus norvegicus',
                  'RN6_ERCC':'Rattus norvegicus',
                           }

    #make list to populate and write to table
    genialisTable=[]
    
    #add the header to the list
    genialisTable.append(genialis_header)
    
    #parse the data table
    dataTable = utils.parseTable(dataFile,'\t')


    #Iterates through the table to select core information for annotation table
    for line in dataTable[1:]:
        name = line[3]
        if paired==False:
            print(line)
            fastq1 = line[8]
            fastq2 = ''
            pair=0
        elif paired==True:
            foo = line[8].split('::')
            print(foo)
            fastq1 = foo[0]
            fastq2 = foo[1]
            pair=1
        if organism=='':
            organism=genomeDict[line[2].upper()]

        new_line = [name,fastq1,fastq2,seqType,pair,collection,annotator,source,organism,strain,tissue,age,
                    genotype,molecule,libraryStrategy,exPro,libraryConst,other1,other2]

        genialisTable.append(new_line)

    #Writes table to output
    utils.unParseTable(genialisTable,outFilePath,'\t')



#==========================================================================
#===================LOADING THE MASTER DATA TABLE==========================
#==========================================================================

def loadDataTable(dataFile):

    if type(dataFile) == str:
        dataTable = parseTable(dataFile,'\t')
    else:
        dataTable = list(dataFile)
    #first check to make sure the table is formatted correctly
    for line in dataTable:
        #print(line)
        if len(line) != 9:
            print('this line did not pass')
            print(line)
            dataTable = formatDataTable(dataFile)
            break
        
    dataDict = {}

    
    
    for line in dataTable[1:]:
        
        dataDict[line[3]] = {}

        dataDict[line[3]]['folder'] = formatFolder(line[0],False)
        dataDict[line[3]]['uniqueID'] = line[1]
        dataDict[line[3]]['genome']=string.upper(line[2])
        genome = line[2]

        dataDict[line[3]]['sam'] = line[0]+line[1]+'.'+genome+'.bwt.sam'
        dataDict[line[3]]['ylf'] = line[0]+line[1]+'.'+genome+'.bwt.ylf'
        dataDict[line[3]]['enriched'] = line[5]
        dataDict[line[3]]['background'] = line[4]
        dataDict[line[3]]['enrichedMacs'] = line[6]
        color_string = string.replace(line[7],'"','')
        dataDict[line[3]]['color'] = color_string
        dataDict[line[3]]['fastq']=line[8]


        #figure out which bam convention we are using
        #default will be new convention
        #look in the bamFolder for all bams that might fit the bill
        bamFolder = "%s" % (line[0])
        bamFileList = [x for x in os.listdir(bamFolder) if len(x) > 0 and x[0] != '.']

        bamFileCandidates = [x for x in bamFileList if x.count(line[1]) == 1 and x.split('.')[-1] =='bam' and x.count('bai') ==0]
        if len(bamFileCandidates) == 0:
            print("UNABLE TO FIND A BAM FILE IN %s WITH UNIQUE ID %s" % (bamFolder,line[1]))
            fullBamPath = ''
        elif len(bamFileCandidates) > 1:
            print("MUTLIPLE BAM FILES IN %s WITH UNIQUE ID %s. NO BAM ASISGNED" % (bamFolder,line[1]))
            print(bamFileCandidates)
            fullBamPath = ''
        else:
            bamFile = bamFileCandidates[0]
            fullBamPath = os.path.abspath('%s%s' % (bamFolder,bamFile))
            fullBaiPath = fullBamPath + '.bai'
            
        if len(fullBamPath) > 0:
            try:
                bam = open(fullBamPath,'r')
                bam.close()
            except IOError:
                print("ERROR: BAM FILE %s DOES NOT EXIST" % (fullBamPath))
                fullBamPath = ''
            try:
                bai = open(fullBaiPath,'r')
                bai.close()
            except IOError:
                print("ERROR: BAM FILE %s DOES NOT HAVE BAI INDEX" % (fullBamPath))
                fullBamPath = ''


        dataDict[line[3]]['bam'] = fullBamPath
            

    return dataDict

def writeDataTable(dataDict,dataFile):
    
    '''
    writes a dataDict to a dataFile
    '''

    newDataTable = [['FILE_PATH','UNIQUE_ID','GENOME','NAME','BACKGROUND','ENRICHED_REGION','ENRICHED_MACS','COLOR','FASTQ_FILE']]

    namesList = dataDict.keys()

    namesList.sort()

    for name in namesList:

        file_path = dataDict[name]['folder']
        uniqueID = dataDict[name]['uniqueID']
        genome = dataDict[name]['genome']
        background = dataDict[name]['background']
        enriched = dataDict[name]['enriched']
        macsEnriched = dataDict[name]['enrichedMacs']
        color = dataDict[name]['color']
        fastq = dataDict[name]['fastq']

        newLine = [file_path,uniqueID,genome,name,background,enriched,macsEnriched,color,fastq]
        newDataTable.append(newLine)

    unParseTable(newDataTable,dataFile,'\t')


def summary(dataFile,outputFile=''):

    '''
    gives a summary of the data and its completeness
    '''

    dataDict = loadDataTable(dataFile)

    dataList = dataDict.keys()
    dataList.sort()
    output= []
    #for each dataset
    isComplete =True
    print('Summarizing data file: %s' %(dataFile))

    for name in dataList:
        print(name)
        #check for the bam
        try:
            bam = open(dataDict[name]['bam'],'r')
            hasBam = True
        except IOError:
            hasBam = False
            isComplete = False

            
        if hasBam == False:
            print('No .bam file for %s' % (name))
            output.append('No .bam file for %s' % (name))

        #if a background is specifed, make sure it exists
        background_name = dataDict[name]['background']
        if len(background_name) > 0 and background_name != 'NONE':
            if background_name in dataDict.keys() == False:
                print('NO BACKGROUND %s DETECED FOR SAMPLE %s' % (background_name,name))

    if outputFile:
        unParseTable(output,outputFile,'')

    if isComplete:
        print('All datasets accounted for in %s\n\n' % (dataFile))
def makeBamTable(dataFile,output):

    '''
    converts a data table into a bam table for jenkins
    schema = [['SOURCE','CELL_TYPE','GENOME','BAMFILE']]
    '''
    
    #sources are manually curated here in this sourceDict

    sourceDict = {'MM1S':'Multiple Myeloma',
                  'KMS11':'Multiple Myeloma',
                  'LY1': 'Diffuse large B-cell lymphoma',
                  'LY3': 'Diffuse large B-cell lymphoma',
                  'LY4': 'Diffuse large B-cell lymphoma',
                  'LY18': 'Diffuse large B-cell lymphoma',
                  'K422': 'Diffuse large B-cell lymphoma',
                  'PFEIFFER':'Diffuse large B-cell lymphoma',
                  'DHL6':'Diffuse large B-cell lymphoma',
                  'HBL1':'Diffuse large B-cell lymphoma',
                  'TOLEDO':'Diffuse large B-cell lymphoma',
                  'P397':'Diffuse large B-cell lymphoma',
                  'P286':'Diffuse large B-cell lymphoma',
                  'P14A':'Lymph node',
                  'P107A':'Tonsil',
                  'P448': 'Diffuse large B-cell lymphoma',
                  'P265': 'Diffuse large B-cell lymphoma',
                  'CD19': 'B-cell',
                  'PROB': 'Pro B-cell',
                  'SKNAS': 'Neuroblastoma',
                  'BE2C': 'Neuroblastoma',
                  'HSC': 'Hematopoeitic stem cell',
                  'EC': 'Human umbilical cord endothelial cell',
                  '3T3L1': 'Fibroblast',
                  'KBM7': 'Haploid Chronic Myeloid Leukemia',
                  'P493-6': 'Burkitt Lymphoma',
                  'NUT797':'Nut midline carcinoma',
                  'MAC': 'Macrophage',
                  'MMP1': 'Primary multiple myeloma',
                  'MEF': 'Mouse embryonic fibroblast',
                  
                  }
    dataDict= loadDataTable(dataFile)
    
    namesList = dataDict.keys()
    namesList.sort()

    #bamTable = [['SOURCE','CELL_TYPE','GENOME','BAMFILE']]
    bamTable =[]
    for name in namesList:

        cellType = name.split('_')[0]
        if sourceDict.has_key(cellType):
            source = sourceDict[cellType]
        else:
            source = 'Unknown'


        genome = dataDict[name]['genome']
        bamFile = dataDict[name]['bam']

        bamTable.append([genome,source,cellType,name,bamFile])


    #sort by source
    sourceOrder = order([x[1] for x in bamTable])
    sortedTable = [['GENOME','SOURCE','CELL_TYPE','NAME','BAMFILE']] + [bamTable[i] for i in sourceOrder]

    unParseTable(sortedTable,output,'\t')
    



#mapping the data

#-------------------------------------------------------------------------#
#                                                                         #
#                          TONY DATABASE TOOLS                            #
#                                                                         #
#-------------------------------------------------------------------------#

#==========================================================================
#=======================INTERACTING WITH TONY==============================
#==========================================================================

def getTONYInfo(uniqueID,column =''):

    '''
    returns a TONY db column parameter
    '''
    column = str(column)

    if len(column) == 0:
        cmd = 'perl /ark/tony/admin/getDB_Data.pl -h'
        os.system(cmd)
        return
    else:
        cmd = 'perl /ark/tony/admin/getDB_Data.pl -i %s -c %s -o TAB' % (uniqueID,column)
        output = subprocess.Popen(cmd,stdin = subprocess.PIPE,stderr = subprocess.PIPE,stdout = subprocess.PIPE,shell = True)

        outputLines = output.stdout.readlines()
        output.stdout.close()
        if outputLines.count('not a valid ID') == 1:
            print("NO RECORD OF %s" % (uniqueID))
            return False
        else:
            return outputLines[1].rstrip().split('\t')[-1]
    
#-------------------------------------------------------------------------#
#                                                                         #
#                                ALIGNMENT                                #
#                                                                         #
#-------------------------------------------------------------------------#


#==========================================================================
#===================CALLING BOWTIE TO MAP DATA=============================
#==========================================================================

def makeBowtieBashJobs(dataFile,namesList = [],launch=True,overwrite=False,paramString=''):

    '''
    makes a mapping bash script and launches 
    '''

    #hardCoded index locations
    dataDict = loadDataTable(dataFile)

    #print(dataDict)
    if len(namesList) == 0:
        namesList = dataDict.keys()
    namesList.sort()
    
    
    for name in namesList:
    
        fastqFile = dataDict[name]['fastq']
        #paired end files will be comma separated
        if fastqFile.count(fastqDelimiter) == 1:
            pairedEnd = True
        elif fastqFile.count(fastqDelimiter) > 1:
            print("UNABLE TO PARSE OUT FASTQ FILES FOR %s" % (name))
        else:
            pairedEnd = False
        genome = dataDict[name]['genome']
        
        #get the unique ID
        uniqueID = dataDict[name]['uniqueID']

        #see if the dataset is already entered into TONY
        #get the parent tony folder
        #tonyFolder = getTONYInfo(uniqueID,column = 30)
        # print(tonyFolder)
        # if tonyFolder:
        #     outputFolder = tonyFolder
        # else:
        #     outputFolder = dataDict[name]['folder']

        outputFolder = dataDict[name]['folder']
        outputFolder = formatFolder(outputFolder,create=True)

        #setting up the folder for linking
        linkFolder = '/grail/bam/%s/' % (string.lower(genome))

        #decide whether or not to run
        try:
            foo = open(dataDict[name]['bam'],'r')
            if not overwrite:
                print('BAM file already exists for %s. OVERWRITE = FALSE' % (name))
                sys.exit()
            else:
                run = True
        except IOError:
            print('no bam file found for %s, making mapping bash script' % (name))
            run = True

        if run:
    
            cmd = "python %s/callBowtie2.py -f %s -g %s -u %s -o %s --link-folder %s" % (whereAmI,fastqFile,genome,uniqueID,outputFolder,linkFolder)
            
            #add the param string
            cmd += " --param '%s'" % (paramString)
            
            if pairedEnd:
                cmd += ' -p'
                
            print(cmd)
            os.system(cmd)
            if launch:
                time.sleep(1)
                cmd = "bash %s%s_bwt2.sh &" % (outputFolder,uniqueID)
                os.system(cmd)

def makeBowtieBashJobsSlurm(dataFile,namesList = [],launch=True,overwrite=False,pCount=1,paramString=' '):

    '''
    makes a mapping bash script and launches
    '''
    if paramString.count('-p') >0:
        print('Error: specify processor count in the pCount argument')
        sys.exit()
    if paramString[0] == ' ':
        paramString = paramString[1:] 
    paramString += ' -p %s' % (pCount)
    print('Using a parameter string of:')
    print(paramString)
    #hardCoded index locations
    dataDict = loadDataTable(dataFile)

    #print(dataDict)
    if len(namesList) == 0:
        namesList = dataDict.keys()
    namesList.sort()


    for name in namesList:

        fastqFile = dataDict[name]['fastq']
        #paired end files will be comma separated
        if fastqFile.count(fastqDelimiter) == 1:
            pairedEnd = True
        elif fastqFile.count(fastqDelimiter) > 1:
            print("UNABLE TO PARSE OUT FASTQ FILES FOR %s" % (name))
        else:
            pairedEnd = False
        genome = dataDict[name]['genome']

        #get the unique ID
        uniqueID = dataDict[name]['uniqueID']

        outputFolder = dataDict[name]['folder']
        outputFolder = formatFolder(outputFolder,create=True)

        #setting up the folder for linking
        linkFolder = '/storage/cylin/grail/bam/%s/' % (string.lower(genome))

        #decide whether or not to run
        try:
            foo = open(dataDict[name]['bam'],'r')
            if not overwrite:
                print('BAM file already exists for %s. OVERWRITE = FALSE' % (name))
                run = False
                pass
            else:
                run = True
        except IOError:
            print('no bam file found for %s, making mapping bash script' % (name))
            run = True

        if run:

            cmd = "python2 %s/callBowtie2Slurm.py -f %s -g %s -u %s -o %s --link-folder %s" % (whereAmI,fastqFile,genome,uniqueID,outputFolder,linkFolder)

            #add the param string
            cmd += " --param '%s'" % (paramString)

            if pairedEnd:
                cmd += ' -p'

            print(cmd)
            os.system(cmd)
            if launch:
                time.sleep(1)
                cmd = "sbatch -n %s --mem 32768 %s%s_bwt2.sh &" % ((pCount+2),outputFolder,uniqueID)
                os.system(cmd)
                #change it to sbatch and off you go


#================================================================
#===================SPLITTING CHIPRX BAMS========================
#================================================================

def splitChipRXBams(dataFile,genome1='',genome2='',namesList=[],header1='',header2=''):
    #genome1 is primary genome
    #genome2 is secondary genome
    #header 1 is the path to the header for genome1 sam files
    #header 2 is the path to the header for genome2 sam files



    dataDict = loadDataTable(dataFile)

    if len(namesList) == 0:
        namesList = dataDict.keys()
    namesList.sort()

    for name in namesList: 
        bam_dir = dataDict[name]['folder']
        uniqueID=dataDict[name]['uniqueID']


        #open the bashfile to write to
        bashFileName = "%s%s_splitBamSlurm.sh" % (bam_dir,uniqueID)
        bashFile = open(bashFileName,'w')

        #shebang
        bashFile.write('#!/usr/bin/bash\n')

        #sbatch funky junk
        cmd = '#SBATCH --output=/storage/cylin/grail/slurm_out/splitChipRXBams_%s' % (uniqueID) + '_%j.out # Standard output and error log'
        bashFile.write(cmd+'\n')
        cmd = '#SBATCH -e /storage/cylin/grail/slurm_out/splitChipRXBams_%s' % (uniqueID) + '_%j.err # Standard output and error log'
        bashFile.write(cmd+'\n')
        bashFile.write('\n\n\n')

        

        cmd = 'samtools view -h -o %s%s.sam %s%s.%s_%s.bwt2.sorted.bam' % (bam_dir,uniqueID,bam_dir,uniqueID,genome1,genome2)
        bashFile.write(cmd+'\n')

        outSam = '%s%s.sam' % (bam_dir,uniqueID)
        split_cmd_1 = 'samtools view -S -F 4 %s | awk -vOFS="\t" '% (outSam)
        split_cmd_2 = '{sindex=index($3,"_%s"); if(sindex >0) {$3=substr($3,1,sindex-1);  print $0 > "%s%s.%s.sam";} else { print $0 > "%s%s.%s.sam";}}' % (genome2,bam_dir,uniqueID,genome2,bam_dir,uniqueID,genome1)
        cmd= split_cmd_1 + ' ' +'\'%s\'' % (split_cmd_2)
        bashFile.write(cmd+'\n')

        genome1_sam = '%s%s.%s.sam' % (bam_dir,uniqueID,genome1)
        genome2_sam = '%s%s.%s.sam' % (bam_dir,uniqueID,genome2)
        outbam1 = '%s%s.%s.bam' % (bam_dir,uniqueID,genome1)
        outbam2 = '%s%s.%s.bam' % (bam_dir,uniqueID,genome2)
        cmd1 = 'cat %s %s | samtools view -bS - > %s' % (header1,genome1_sam, outbam1)
        cmd2 = 'cat %s %s | samtools view -bS - > %s' % (header2,genome2_sam, outbam2)
        bashFile.write(cmd1 + '\n' + cmd2 + '\n\n')

        cmd = 'samtools sort \'%s\' \'%s%s.%s.sorted\''  % (outbam1,bam_dir,uniqueID,genome1)
        bashFile.write(cmd + '\n')
        cmd = 'samtools sort \'%s\' \'%s%s.%s.sorted\''  % (outbam2,bam_dir,uniqueID,genome2)
        bashFile.write(cmd + '\n')
        cmd = 'samtools index %s%s.%s.sorted.bam' % (bam_dir,uniqueID,genome1)
        bashFile.write(cmd + '\n')
        cmd = 'samtools index %s%s.%s.sorted.bam' % (bam_dir,uniqueID,genome2)
        bashFile.write(cmd + '\n')

        bashFile.close()


#===========================================================================
#=========================ChIP RX SCALE FACTORS=============================
#===========================================================================




def writeScaleFactors(dataFile,namesList=[],output='',genome1='',genome2=''):

    '''
    creates a table of scale factors based on the rx genome read depth
    '''

    #first set up the output folder
    #rpm scale factor is what the rpm/bp should be MULTIPLIED by
    #mouse mapped reads give the denominator for what raw r/bp should be divided by
    genome1_map_string = '%s_MAPPED_READS' % (genome1)
    genome2_map_string = '%s_MAPPED_READS' % (genome2)   
    outputTable = [['NAME',genome1_map_string,genome2_map_string,'RPM_SCALE_FACTOR']]


    dataDict=loadDataTable(dataFile)
    if len(namesList) == 0:
        namesList = [name for name in dataDict.keys()]
    namesList.sort()
    print('scaling the following datasets')


    for name in namesList:

        print('WORKING ON %s' % (name))
        bam_path = dataDict[name]['bam']
        bam = utils.Bam(bam_path)
        bam_mmr = float(bam.getTotalReads())/1000000
        scale_path = string.replace(bam_path,genome1,genome2)
        scaleBam = utils.Bam(scale_path)
        scale_mmr = float(scaleBam.getTotalReads())/1000000
        print(bam_path)
        print(bam_mmr)
        print(scale_path)
        print(scale_mmr)
        rpm_scale = bam_mmr/scale_mmr
        scale_line = [bam_mmr,scale_mmr,rpm_scale]
        scale_line = [round(x,4) for x in scale_line]
        outputTable.append([name] + scale_line)

    if len(output) == 0:
        return outputTable
    else:
        utils.unParseTable(outputTable,output,'\t')







###
### DEPRECATED FUNCTION FROM TAK
###
# def callBowtie(dataFile,dataList = [],overwrite = False,pairedEnd = False):

#     '''
#     calls bowtie for the dataset names specified. if blank, calls everything
#     '''
    
#     dataDict = loadDataTable(dataFile)

#     if len(dataList) == 0:
#         dataList = dataDict.keys()
    
#     for name in dataList:
        
#         #make sure the output folder exists
#         try:
#             foo = os.listdir(dataDict[name]['folder'])
#         except OSError:
#             print('no output folder %s for dataset %s. Creating directory %s now.' % (dataDict[name]['folder'],name,dataDict[name]['folder']))
#             os.system('mkdir %s' % (dataDict[name]['folder']))

            
#         if string.upper(dataDict[name]['genome']) == 'HG18' or string.upper(dataDict[name]['genome']) == 'MM8':
#             cmd = 'perl /nfs/young_ata/scripts/generateSAM_file.pl -R -G -c 4 -B -i %s -o %s' % (dataDict[name]['uniqueID'],dataDict[name]['folder'])
#         elif string.upper(dataDict[name]['genome']) == 'RN5':
#             cmd = 'perl /nfs/young_ata/CYL_code/generateSam_rat.pl -R -c 4 -B -i %s -o %s' % (dataDict[name]['uniqueID'],dataDict[name]['folder'])
#         else:
#             cmd = 'perl /nfs/young_ata/scripts/generateSAM_file.pl -R -c 4 -B -i %s -o %s' % (dataDict[name]['uniqueID'],dataDict[name]['folder'])
        
#         #first check if the bam exists
#         if overwrite:
#             print(cmd)
#             os.system(cmd)
#         else:
#             try:
#                 foo = open(dataDict[name]['bam'],'r')
#                 print('BAM file already exists for %s. OVERWRITE = FALSE' % (name))
#             except IOError:
#                 print('no bam file found for %s, mapping now' % (name))
#                 print(cmd)
#                 os.system(cmd)



#==========================================================================
#===================GETTING MAPPING STATS==================================
#==========================================================================

def bowtieStats(dataFile,namesList=[]):

    '''
    gets stats from the bowtie out file for each bam
    does not assume bowtie output format is always the same
    '''
    dataDict = loadDataTable(dataFile)
    print('DATASET\t\tSEED LENGTH\tTOTAL READS (MILLIONS)\tALIGNED\tFAILED\tMULTI')
    if len(namesList) == 0:
        namesList= dataDict.keys()
    for name in namesList:

        readLength,readCount,alignedReads,failedReads,multiReads = False,False,False,False,False
        bowtieOutFile = dataDict[name]['folder'] + dataDict[name]['uniqueID']+ '.bowtie.output'
        try:
            bowtieOut = open(bowtieOutFile,'r')
        except IOError:
            print('NO BOWTIE OUTPUT FOR %s' % (name))
            continue
        for line in bowtieOut:
            #get the read length
            if line[0] == '':
                continue
            if line[0:4] == 'File':
                readLength = line.split('=')[1][1:-1]

            if line[0] == '#' and line.count('reads') == 1:
                
                if line.count('processed') == 1:

                    readCount = line.split(':')[1][1:-1]
                    readCount = round(float(readCount)/1000000,2)
                if line.count('reported alignment') == 1:
                    alignedReads = line[-9:-1]

                if line.count('failed') == 1:
                    failedReads = line[-9:-1]
                if line.count('suppressed') == 1:
                    multiReads = line[-9:-1]

        if readLength and readCount and alignedReads and failedReads and multiReads:
            print('%s\t\t%s\t%s\t%s\t%s\t%s' % (name,readLength,readCount,alignedReads,failedReads,multiReads))
        else:
            print('%s\tNO DATA AVAILABLE' % (name))


        
        

#==========================================================================
#===================MERGE BAMS=============================================
#==========================================================================

def mergeBams(dataFile,mergedName,namesList,color='',background =''):

    '''
    merges a set of bams and adds an entry to the dataFile
    '''
    dataDict = loadDataTable(dataFile)

    #make an .sh file to call the merging
    
    bamFolder = dataDict[namesList[0]]['folder']
    genome = string.lower(dataDict[namesList[0]]['genome'])
    if color == '':
        color = dataDict[namesList[0]]['color']

    if background == '':
        background = dataDict[namesList[0]]['background']

        
    
    if bamFolder[-1] != '/':
        bamFolder+='/'
    mergedFullPath = bamFolder+mergedName +'.'+genome+'.bwt.sorted.bam'

    mergeShFile = open(mergedFullPath+'.sh','w')

    mergeShFile.write('cd %s\n' % (bamFolder))

    cmd1 = '%s merge %s ' % (samtoolsString,mergedFullPath)

    for name in namesList:

        cmd1+= ' %s' % (dataDict[name]['bam'])

    mergeShFile.write(cmd1+'\n')

    cmd2 = '%s sort %s %s' % (samtoolsString,mergedFullPath,bamFolder+mergedName+'.'+genome+'.bwt.sorted')
    mergeShFile.write(cmd2+'\n')

    cmd3 = '%s index %s' % (samtoolsString,mergedFullPath)
    mergeShFile.write(cmd3+'\n')

    mergeShFile.close()

    runCmd = " 'bash %s'" % (mergedFullPath+'.sh')
    os.system(runCmd)

    dataTable = parseTable(dataFile,'\t')

    dataTable.append([bamFolder,mergedName,genome,mergedName,background,'','',color])
    unParseTable(dataTable,dataFile,'\t')
    formatDataTable(dataFile)


#==========================================================================
#===============================FILTERING BAMS=============================
#==========================================================================


def filterBams(dataFile,namesList = [],tempFolder = '/raider/BOWTIE_TEMP/',bamFolder='/grail/bam/filtered/'):

    '''
    for all datasets, will run the bam filtering, drop the filtered bam into the correct TONY folder
    and symlink into the bamFolder
    must run as admin to write to folder
    '''

    #check that the temp folder and bam folder exist
    if not formatFolder(tempFolder,False):
        print("ERROR: UNABLE TO FIND TEMPORARY FOLDER %s" % (tempFolder))
        sys.exit()
    else:
        tempFolder = formatFolder(tempFolder,False)
    if not formatFolder(bamFolder,False):
        print("ERROR: UNABLE TO FIND BAM FOLDER %s" % (bamFolder))
        sys.exit()
    else:
        bamFolder = formatFolder(bamFolder,False)
    
    dataDict = loadDataTable(dataFile)

    if len(namesList) == 0:
        namesList = dataDict.keys()

    
    for name in namesList:

        uniqueID = dataDict[name]['uniqueID']
        genome = string.lower(dataDict[name]['genome'])
        
        #check to see if this has already been done
        if dataDict[name]['folder'] == "%s%s/" % (bamFolder,genome):
            continue

        #change the bam Folder to the location of the filtered bam
        dataDict[name]['folder'] = "%s%s/" % (bamFolder,genome)
        #get the original bam
        bamFile = getTONYInfo(uniqueID,column =47)

        #get the parent folder
        parentFolder = getParentFolder(bamFile)

        #make a temp folder
        tempOutFolder = tempFolder + uniqueID +'/'
        tempOutFolder = formatFolder(tempOutFolder,True)

        #set up a bashFile
        bashFileName = '%s%s_filterBam.sh' % (tempOutFolder,uniqueID)
        bashFile = open(bashFileName,'w')
        bashFile.write('#!usr/bin/bash\n')
        
        #cd into temp directory
        bashFile.write('cd %s\n' % (tempOutFolder))
        
        #set up the filter command
        filterCmd = 'python /ark/home/cl512/src/scripts-for-project/2-mito-filter-sort-dedupe.py -t 4 -i %s -o %s' % (bamFile,tempOutFolder)
        bashFile.write(filterCmd + '\n')

        #set up the command to move and then sym link
        #mv command
        mvBamCommand = 'mv *noChrM.sort.rmdup.bam %s%s_%s.noChrM.fix.rmdup.sorted.bam' % (parentFolder,uniqueID,genome)
        bashFile.write(mvBamCommand + '\n')

        mvBaiCommand = 'mv *noChrM.sort.rmdup.bam.bai %s%s_%s.noChrM.fix.rmdup.sorted.bam.bai' % (parentFolder,uniqueID,genome)
        bashFile.write(mvBaiCommand + '\n')

        #now symlink
        #check to make sure the genome folder exist
        if not formatFolder(bamFolder+genome,False):
            print("NO FOLDER EXISTS IN %s FOR BUILD %s" % (bamFolder,genome))
            print("ATTEMPTING TO CREATE FOLDER %s%s/" % (bamFolder,genome))
        symFolder = formatFolder(bamFolder+genome,True)

        symBamCommand = 'ln -s %s%s_%s.noChrM.fix.rmdup.sorted.bam %s' % (parentFolder,uniqueID,genome,symFolder)
        bashFile.write(symBamCommand + '\n')

        symBaiCommand = 'ln -s %s%s_%s.noChrM.fix.rmdup.sorted.bam.bai %s' % (parentFolder,uniqueID,genome,symFolder)
        bashFile.write(symBaiCommand + '\n')

        #now the cleanup
        cleanCommand = 'rm -rf %s' % (tempOutFolder)
        bashFile.write(cleanCommand + '\n')
        bashFile.close()
        runCommand = 'bash %s &' %(bashFileName) 
        print("Run command: %s" % (runCommand))
        os.system(runCommand)

    
    writeDataTable(dataDict,dataFile)


#-------------------------------------------------------------------------#
#                                                                         #
#                              PEAK FINDING                               #
#                                                                         #
#-------------------------------------------------------------------------#


#==========================================================================
#===================CALLING MACS ==========================================
#==========================================================================

def callMacsQsub(dataFile,macsFolder,namesList = [],overwrite=False,pvalue='1e-9'):

    '''
    calls the macs error model
    '''
    dataDict = loadDataTable(dataFile)

    if macsFolder[-1] != '/':
        macsFolder+='/'
    formatFolder(macsFolder,True)

    #if no names are given, process all the data
    if len(namesList) == 0:

        namesList = dataDict.keys()

    #a timestamp to name this pipeline batch of files
    timeStamp =datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    #a random integer ticker to help name files
    randTicker = random.randint(0,10000)

    for name in namesList:
        
        #skip if a background set
        if string.upper(dataDict[name]['background']) == 'NONE':
            continue

        #check for the bam
        try:
            bam = open(dataDict[name]['bam'],'r')
            hasBam = True
        except IOError:
            hasBam = False

        #check for background
        try:
            backgroundName = dataDict[name]['background']
            backbroundBam = open(dataDict[backgroundName]['bam'],'r')
            hasBackground = True
        except IOError:
            hasBackground = False

        if not hasBam:
            print('no bam found for %s. macs not called' % (name))
            continue

        if not hasBackground:
            print('no background bam %s found for dataset %s. macs not called' % (backgroundName,name))
            continue
        #make a new folder for every dataset
        outdir = macsFolder+name
        #print(outdir)
        try:
            foo = os.listdir(outdir)
            if not overwrite:
                print('MACS output already exists for %s. OVERWRITE = FALSE' % (name))
                continue
            
        except OSError:
            os.system('mkdir %s' % (outdir))
            
        os.chdir(outdir)
        genome = dataDict[name]['genome']
        #print('USING %s FOR THE GENOME' % genome)
        
        bamFile = dataDict[name]['bam']
        backgroundName =  dataDict[name]['background']
        backgroundBamFile = dataDict[backgroundName]['bam']
        macsString = '/usr/local/python-2.7.2/bin/macs14'
        if string.upper(genome[0:2]) == 'HG':
            cmd = "%s -t %s -c %s -f BAM -g hs -n %s -p %s -w -S --space=50" % (macsString,bamFile,backgroundBamFile,name,pvalue)
        elif string.upper(genome[0:2]) == 'MM':
            cmd = "%s -t %s -c %s -f BAM -g mm -n %s -p %s -w -S --space=50" % (macsString,bamFile,backgroundBamFile,name,pvalue)
        elif string.upper(genome[0:2]) == 'RN':
            cmd = "%s -t %s -c %s -f BAM -g 2000000000 -n %s -p %s -w -S --space=50" % (macsString,bamFile,backgroundBamFile,name,pvalue)
        print(cmd)

        bashFileName = '/mnt/d0-0/share/bradnerlab/src/cl512/temp/macs_%s_%s.sh' % (timeStamp,randTicker)
        bashFile = open(bashFileName,'w')
        bashFile.write("cd %s\n" % (outdir))
        bashFile.write(cmd)
        bashFile.close()
        bashCommand = 'qsub %s' % (bashFileName) 
        print(bashCommand)
        os.system(bashCommand)

        randTicker+=1



def callMacs(dataFile,macsFolder,namesList = [],overwrite=False,pvalue='1e-9',useBackground =True):

    '''
    calls the macs error model
    '''
    dataDict = loadDataTable(dataFile)

    if macsFolder[-1] != '/':
        macsFolder+='/'
    formatFolder(macsFolder,True)

    #if no names are given, process all the data
    if len(namesList) == 0:

        namesList = dataDict.keys()


    for name in namesList:
        
        #skip if a background set
        if useBackground and string.upper(dataDict[name]['background']) == 'NONE':
            continue

        #check for the bam
        try:
            bam = open(dataDict[name]['bam'],'r')
            hasBam = True
        except IOError:
            hasBam = False

        #check for background
        try:
            backgroundName = dataDict[name]['background']
            backbroundBam = open(dataDict[backgroundName]['bam'],'r')
            hasBackground = True
        except (IOError, KeyError) as e:
            hasBackground = False

        if not hasBam:
            print('no bam found for %s. macs not called' % (name))
            continue

        if useBackground and hasBackground == False:
            print('no background bam %s found for dataset %s. macs not called' % (backgroundName,name))
            continue
        #make a new folder for every dataset
        outdir = macsFolder+name
        #print(outdir)
        try:
            foo = os.listdir(outdir)
            if not overwrite:
                print('MACS output already exists for %s. OVERWRITE = FALSE' % (name))
                continue
            
        except OSError:
            
            os.system('mkdir %s' % (outdir))
            
        os.chdir(outdir)
        genome = dataDict[name]['genome']
        #print('USING %s FOR THE GENOME' % genome)
        if useBackground == True:
            bamFile = dataDict[name]['bam']
            backgroundName =  dataDict[name]['background']
            backgroundBamFile = dataDict[backgroundName]['bam']
            if string.upper(genome[0:2]) == 'HG':
                cmd = "macs14 -t %s -c %s -f BAM -g hs -n %s -p %s -w -S --space=50 &" % (bamFile,backgroundBamFile,name,pvalue)
            elif string.upper(genome[0:2]) == 'MM':
                cmd = "macs14 -t %s -c %s -f BAM -g mm -n %s -p %s -w -S --space=50 &" % (bamFile,backgroundBamFile,name,pvalue)
            elif string.upper(genome[0:2]) == 'RN':
                cmd = "macs14 -t %s -c %s -f BAM -g 2000000000 -n %s -p %s -w -S --space=50 &" % (bamFile,backgroundBamFile,name,pvalue)

            elif string.upper(genome[0:2]) == 'DA':
                cmd = "macs14 -t %s -c %s -f BAM -g 1000000000 -n %s -p %s -w -S --space=50 &" % (bamFile,backgroundBamFile,name,pvalue)

        if useBackground == False:
            bamFile = dataDict[name]['bam']


            if string.upper(genome[0:2]) == 'HG':
                cmd = "macs14 -t %s -f BAM -g hs -n %s -p %s -w -S --space=50 &" % (bamFile,name,pvalue)
            elif string.upper(genome[0:2]) == 'MM':
                cmd = "macs14 -t %s -f BAM -g mm -n %s -p %s -w -S --space=50 &" % (bamFile,name,pvalue)
            elif string.upper(genome[0:2]) == 'RN':
                cmd = "macs14 -t %s -f BAM -g 2000000000 -n %s -p %s -w -S --space=50 &" % (bamFile,name,pvalue)
            elif string.upper(genome[0:2]) == 'DA':
                cmd = "macs14 -t %s -f BAM -g 1000000000 -n %s -p %s -w -S --space=50 &" % (bamFile,name,pvalue)



        print(cmd)
        os.system(cmd)



def callMacsSlurm(dataFile,macsFolder,namesList = [],overwrite=False,pvalue='1e-9',useBackground =True,launch=True):

    '''
    calls the macs error model
    '''
    dataDict = loadDataTable(dataFile)

    if macsFolder[-1] != '/':
        macsFolder+='/'
    formatFolder(macsFolder,True)

    #if no names are given, process all the data
    if len(namesList) == 0:

        namesList = dataDict.keys()


    for name in namesList:

        #skip if a background set
        if useBackground and string.upper(dataDict[name]['background']) == 'NONE':
            continue

        #check for the bam
        try:
            bam = open(dataDict[name]['bam'],'r')
            hasBam = True
        except IOError:
            hasBam = False

        #check for background
        try:
            backgroundName = dataDict[name]['background']
            backbroundBam = open(dataDict[backgroundName]['bam'],'r')
            hasBackground = True
        except (IOError, KeyError) as e:
            hasBackground = False

        if not hasBam:
            print('no bam found for %s. macs not called' % (name))
            continue

        if useBackground and hasBackground == False:
            print('no background bam %s found for dataset %s. macs not called' % (backgroundName,name))
            continue
        #make a new folder for every dataset
        outdir = macsFolder+name
        #print(outdir)
        try:
            foo = os.listdir(outdir)
            print('am i checking directory')
            print(foo)
            if not overwrite:
                print('MACS output already exists for %s. OVERWRITE = FALSE' % (name))
                continue

        except OSError:
            print('am i being run')
            os.system('mkdir %s' % (outdir))

        os.chdir(outdir)
        #create a bashfile in the outdir
        #that has the name of the dataset which is a parameter
        bashFileName = '%s/macs14_%s.sh' %(outdir,name)
        bashFile = open(bashFileName,'w')
        #shebang
        bashFile.write('#!/usr/bin/bash\n')
        #write the bash cmd w/ any appropriate slurmy slurm stuff
        #get a time stamp
        ts = time.time()
        timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%Hh%Mm%Ss')
        cmd = '#SBATCH --output=/storage/cylin/grail/slurm_out/macs_%s_%s' % (name,timestamp) + '_%j.out # Standard output and error log'
        bashFile.write(cmd+'\n')
        cmd = '#SBATCH -e /storage/cylin/grail/slurm_out/macs_%s_%s' % (name,timestamp) + '_%j.err # Standard output and error log'
        bashFile.write(cmd+'\n')

        cmd = 'pwd; hostname; date'
        bashFile.write(cmd+'\n\n\n\n')

        genome = dataDict[name]['genome']
        #print('USING %s FOR THE GENOME' % genome)
        if useBackground == True:
            bamFile = dataDict[name]['bam']
            backgroundName =  dataDict[name]['background']
            backgroundBamFile = dataDict[backgroundName]['bam']
            if string.upper(genome[0:2]) == 'HG':
                cmd = "macs14 -t %s -c %s -f BAM -g hs -n %s -p %s -w -S --space=50" % (bamFile,backgroundBamFile,name,pvalue)
            elif string.upper(genome[0:2]) == 'MM':
                cmd = "macs14 -t %s -c %s -f BAM -g mm -n %s -p %s -w -S --space=50" % (bamFile,backgroundBamFile,name,pvalue)
            elif string.upper(genome[0:2]) == 'RN':
                cmd = "macs14 -t %s -c %s -f BAM -g 2000000000 -n %s -p %s -w -S --space=50" % (bamFile,backgroundBamFile,name,pvalue)

            elif string.upper(genome[0:2]) == 'DA':
                cmd = "macs14 -t %s -c %s -f BAM -g 1000000000 -n %s -p %s -w -S --space=50" % (bamFile,backgroundBamFile,name,pvalue)

        if useBackground == False:
            bamFile = dataDict[name]['bam']

            if string.upper(genome[0:2]) == 'HG':
                cmd = "macs14 -t %s -f BAM -g hs -n %s -p %s -w -S --space=50" % (bamFile,name,pvalue)
            elif string.upper(genome[0:2]) == 'MM':
                cmd = "macs14 -t %s -f BAM -g mm -n %s -p %s -w -S --space=50" % (bamFile,name,pvalue)
            elif string.upper(genome[0:2]) == 'RN':
                cmd = "macs14 -t %s -f BAM -g 2000000000 -n %s -p %s -w -S --space=50" % (bamFile,name,pvalue)
            elif string.upper(genome[0:2]) == 'DA':
                cmd = "macs14 -t %s -f BAM -g 1000000000 -n %s -p %s -w -S --space=50" % (bamFile,name,pvalue)

        print(cmd)
        #bashFile.write('sbatch ' + cmd)
        bashFile.write(cmd)
        bashFile.close()
        print('bash file: %s' % (bashFileName))
        #return the path to teh bash file and launch if you want

        if launch == True:
            print('I AM LAUNCHING JOB')
            os.system('sbatch -n 2 --mem 32768 %s' %(bashFileName))
        #return the path to teh bash file and launch if you want



def callMacs2(dataFile,macsFolder,namesList = [],broad=True,noBackground = False,pairedEnd = False,overwrite=False,pvalue='1e-9'):

    '''
    calls the macs error model
    '''
    dataDict = loadDataTable(dataFile)

    if macsFolder[-1] != '/':
        macsFolder+='/'
    formatFolder(macsFolder,True)

    #if no names are given, process all the data
    if len(namesList) == 0:

        namesList = dataDict.keys()


    for name in namesList:
        
        #skip if a background set
        if string.upper(dataDict[name]['background']) == 'NONE' and noBackground == False:
            continue

        #check for the bam
        try:
            bam = open(dataDict[name]['bam'],'r')
            hasBam = True
        except IOError:
            hasBam = False

        #check for background

        if noBackground == False:
            try:
                backgroundName = dataDict[name]['background']
                backbroundBam = open(dataDict[backgroundName]['bam'],'r')
                hasBackground = True
            except IOError:
                hasBackground = False
        else:
            hasBackground = False

        if not hasBam:
            print('no bam found for %s. macs not called' % (name))
            continue

        if not hasBackground and noBackground == False:
            print('no background bam %s found for dataset %s. macs not called' % (backgroundName,name))
            continue
        #make a new folder for every dataset
        outdir = macsFolder+name
        #print(outdir)
        try:
            foo = os.listdir(outdir)
            if not overwrite:
                print('MACS output already exists for %s. OVERWRITE = FALSE' % (name))
                continue
            
        except OSError:
            
            os.system('mkdir %s' % (outdir))
            
        os.chdir(outdir)

        bamFile = dataDict[name]['bam']

        genome = dataDict[name]['genome']
        #print('USING %s FOR THE GENOME' % genome)
        
        #setting parameters
        if string.upper(genome[0:2]) == 'HG':
            genomeString = 'hs'

        if pairedEnd == True:
            fileType = 'BAMPE'
        else:
            fileType = 'BAM'

        if broad == True:
            broadCall = '--broad'
        else:
            broadCall = ''
            
        if noBackground:
            cmd = "macs2 callpeak -t %s -f %s -g %s --outdir %s -n %s -p %s %s &" % (bamFile,fileType,genomeString,outdir,name,pvalue,broadCall)
        else:
            backgroundName =  dataDict[name]['background']
            backgroundBamFile = dataDict[backgroundName]['bam']
            cmd = "macs2 callpeak -t %s -c %s -f %s -g %s --outdir %s -n %s -p %s %s &" % (bamFile,backgroundBamFile,fileType,genomeString,outdir,name,pvalue,broadCall)


        print(cmd)
        os.system(cmd)




#==========================================================================
#===================FORMAT MACS OUTPUT=====================================
#==========================================================================

def formatMacsOutput(dataFile,macsFolder,macsEnrichedFolder,wiggleFolder,wigLink='',useBackground=True):

    dataDict = loadDataTable(dataFile)

    if macsFolder[-1] != '/':
        macsFolder+='/'

    if macsEnrichedFolder[-1] != '/':
        macsEnrichedFolder+='/'

    #make an output directory for the macsEnriched
    formatFolder(macsEnrichedFolder,True)
    #make an output directory for the wiggles
    formatFolder(wiggleFolder,True)


        
    
    namesList = dataDict.keys()
    namesList.sort()
    dataTable = parseTable(dataFile,'\t')

    if len(wigLink) > 0:
        if wigLink[-1] != '/':
            wigLink+='/'
        genome = string.lower(dataDict[namesList[0]]['genome'])
        cellType = string.lower(namesList[0].split('_')[0])
        wigLinkFolder = '%s%s/%s' % (wigLink,genome,cellType)
        formatFolder(wigLinkFolder,True)


    newDataTable = [dataTable[0]]

    dataTableRows = [line[3] for line in dataTable]
    for name in namesList:

        if len(wigLink) > 0:
            if wigLink[-1] != '/':
                wigLink+='/'
            genome = string.lower(dataDict[namesList[0]]['genome'])
            cellType = string.lower(namesList[0].split('_')[0])
            wigLinkFolder = '%s%s/%s' % (wigLink,genome,cellType)
            formatFolder(wigLinkFolder,True)


        outwiggleFileName = '%s%s_treat_afterfiting_all.wig.gz' % (wiggleFolder,name)
        #find the row in the old dataTable
        print('looking for macs output for %s' % (name))
        i = dataTableRows.index(name)
        #skip if a background set
        if useBackground and string.upper(dataDict[name]['background']) == 'NONE':
            newLine = list(dataTable[i])
            newLine[6] = 'NONE'
            newDataTable.append(newLine)
            continue
        #check to see if the wiggle exists in the correct folder
        try:
            outwiggle = open('%s%s_treat_afterfiting_all.wig.gz' % (wiggleFolder,name),'r')
            outwiggle.close()
        except IOError:
            #move the wiggle and add the itemrgb line
            if dataDict[name]['color'] == '0,0,0':
                print('moving wiggle for %s over without changing color' % (name))
                cmd = 'mv %s%s/%s_MACS_wiggle/treat/%s_treat_afterfiting_all.wig.gz %s' % (macsFolder,name,name,name,outwiggleFileName)
                os.system(cmd)
            else:
                try:
                    #print(name)
                    print('for dataset %s going from %s%s/%s_MACS_wiggle/treat/%s_treat_afterfiting_all.wig.gz' % (name,macsFolder,name,name,name))
                    wiggle = open('%s%s/%s_MACS_wiggle/treat/%s_treat_afterfiting_all.wig.gz' % (macsFolder,name,name,name),'r')
                    #print(name)
                    print('and writing to %s%s_treat_afterfiting_all.wig.gz' % (wiggleFolder,name))
                
                    outwiggle = open(outwiggleFileName,'w')

                    print('writing new wiggle with color line for %s' % (name))
                    header = wiggle.readline().rstrip()
                    color = dataDict[name]['color']
                    header = header + ' itemRgb="On" color="%s"' % (color) + '\n'
                    outwiggle.write(header)
                    for line in wiggle:
                        outwiggle.write(line)
                    outwiggle.close()
                    wiggle.close()

                
                except IOError:
                    print('WARNING: NO MACS WIGGLE FOR %s' %(name)) 

        
        if len(wigLink) > 0:
            #time.sleep(10)
            print('Creating symlink for dataset %s in %s' % (name,wigLinkFolder))
            os.system('cd %s' % (wigLinkFolder))
            os.chdir(wigLinkFolder)
            print('cd %s' % (wigLinkFolder))
            wigLinkFileName = '%s_treat_afterfitting_all.wig.gz' % (name)
            print('ln -s %s %s' % (outwiggleFileName,wigLinkFileName))
            os.system('ln -s %s %s' % (outwiggleFileName,wigLinkFileName))


        #first check if the thing exists
        try:
            foo = open('%s%s_peaks.bed' % (macsEnrichedFolder,name),'r')
            newLine = list(dataTable[i])
            newLine[6] = name+'_peaks.bed'
            newDataTable.append(newLine)
        except IOError:
            #move the bedFile of the peaks
        
            try:
                foo = open('%s%s/%s_peaks.bed' % (macsFolder,name,name),'r')
                cmd = 'mv %s%s/%s_peaks.bed %s' % (macsFolder,name,name,macsEnrichedFolder)
                newLine = list(dataTable[i])
                newLine[6] = name+'_peaks.bed'
                newDataTable.append(newLine)
                print(cmd)
                os.system(cmd)
            except IOError:
                print('WARNING: NO MACS OUTPUT FOR %s' %(name)) 
                newLine = list(dataTable[i])
                newLine[6] = 'NONE'
                newDataTable.append(newLine)


        

    unParseTable(newDataTable,dataFile,'\t')




#==========================================================================
#=================OVERALL WRAPPER FOR RUNNING MACS=========================
#==========================================================================




def run_macs(data_file,projectFolder,macsFolder,macsEnrichedFolder,wiggleFolder,useBackground=True):
    dataDict = loadDataTable(data_file)
    namesList = [name for name in dataDict.keys() if name.upper().count('WCE') ==0 and name.upper().count('INPUT') == 0]
    namesList.sort()
    print(namesList)
    callMacs(data_file,macsFolder,namesList,False,'1e-9',useBackground)
    os.chdir(projectFolder) # the silly call macs script has to change into the output dir
    #so this takes us back to the project folder

    #to check for completeness, we will try to find all of the peak files
    peak_calling_done = False
    while not peak_calling_done:
        dataDict = loadDataTable(data_file)
        namesList = [name for name in dataDict.keys() if name.upper().count('WCE') ==0 and name.upper().count('INPUT') == 0]
        for name in namesList:
            peak_path = '%s%s/%s_summits.bed' % (macsFolder,name,name)
            print('searching for %s' % (peak_path))
            if utils.checkOutput(peak_path,1,180):
                peak_calling_done =True
                print('found %s' % (peak_path))
                continue
            else:
                print('Error: peak calling timed out')
                sys.exit()
    
    #now format the macs output
    print('formatting macs output')
    dataDict = loadDataTable(data_file)
    namesList = [name for name in dataDict.keys() if name.upper().count('WCE') ==0 and name.upper().count('INPUT') == 0]
    formatMacsOutput(data_file,macsFolder,macsEnrichedFolder,wiggleFolder,'',useBackground)
    print('Finished running Macs 1.4.2')




#-------------------------------------------------------------------------#
#                                                                         #
#                              GFF TOOLS                                  #
#                                                                         #
#-------------------------------------------------------------------------#



#==========================================================================
#===================MAKING GFFS OF TSS REGIONS=============================
#==========================================================================
        
def makeGeneGFFs(annotFile,gffFolder,species='HG18'):

    '''
    makes a tss gff with the given window size for all genes in the annotation
    tss +/-5kb
    tss +/-300
    body +300/+3000
    can work on any genome build given the right annot file
    '''


    if gffFolder[-1] != '/':
        gffFolder+='/'

    try:
        foo = os.listdir(gffFolder)
        print('Directory %s already exists. Using %s to store gff' % (gffFolder,gffFolder))
    except OSError:
        cmd = 'mkdir %s' % (gffFolder)
        print('making directory %s to store gffs' % (gffFolder))
        os.system(cmd)


    startDict = makeStartDict(annotFile)

    geneList = startDict.keys()
    print("USING %s genes" % (len(geneList)))
    
    tssLoci = []
    for gene in geneList:
        tssLocus = Locus(startDict[gene]['chr'],startDict[gene]['start'][0]-5000,startDict[gene]['start'][0]+5000,startDict[gene]['sense'],gene)
        tssLoci.append(tssLocus)
    #print(len(tssLoci))
    tssCollection = LocusCollection(tssLoci,500)
    #print(len(tssCollection.getLoci()))
    tssGFF_5kb = locusCollectionToGFF(tssCollection)
    #print(len(tssGFF_5kb))
    tssLoci= [] 
    for gene in geneList:
        tssLocus = Locus(startDict[gene]['chr'],startDict[gene]['start'][0]-1000,startDict[gene]['start'][0]+1000,startDict[gene]['sense'],gene)
        tssLoci.append(tssLocus)
    tssCollection = LocusCollection(tssLoci,500)
    tssGFF_1kb = locusCollectionToGFF(tssCollection)
    
    tssGFF_300 = []
    txnGFF = []

    for line in tssGFF_5kb:
        gene = line[1]
        chrom = startDict[gene]['chr']
        start = startDict[gene]['start'][0]
        end = startDict[gene]['end'][0]
        sense = startDict[gene]['sense']
        name = startDict[gene]['name']
        tssLine = [chrom,gene,'',start-300,start+300,'',sense,'',gene]
        if sense == '+':
            txnLine = [chrom,gene,'',start+300,end+3000,'',sense,'',gene]
        else:
            txnLine = [chrom,gene,'',end-3000,start-300,'',sense,'',gene]

        tssGFF_300.append(tssLine)
        txnGFF.append(txnLine)

     
    unParseTable(tssGFF_5kb,gffFolder + '%s_TSS_ALL_-5000_+5000.gff' % (species),'\t')
    unParseTable(tssGFF_1kb,gffFolder + '%s_TSS_ALL_-1000_+1000.gff' % (species),'\t')
    unParseTable(tssGFF_300,gffFolder + '%s_TSS_ALL_-300_+300.gff' % (species),'\t')
    unParseTable(txnGFF,gffFolder + '%s_BODY_ALL_+300_+3000.gff' % (species),'\t')

#==========================================================================
#===================MAKING GFFS OF CHROMS==================================
#==========================================================================
        
def makeChromGFFs(chromLengthFile,gffFolder,chromList = [],genome='HG18',binSize = 100000,singleGFF = True):

    '''
    makes GFFs of chromosomes, each chrom gets its own gff
    '''
    formatFolder(gffFolder,True)
    chromLengthDict = {}
    chromLengthTable = parseTable(chromLengthFile,'\t')


    genomesList = uniquify([line[2] for line in chromLengthTable])

    for x in genomesList:
        chromLengthDict[string.upper(x)] = {}

    for line in chromLengthTable:
        chromLengthDict[string.upper(line[2])][line[0]] = int(line[4])
    if len(chromList) ==0:
        chromList = chromLengthDict[string.upper(genome)].keys()
        chromList.sort()
    masterGFF = []
    for chrom in chromList:
        if chromLengthDict[string.upper(genome)].has_key(chrom):
            chromGFF = []
            ticker = 1
            for i in range(1,chromLengthDict[genome][chrom],int(binSize)):
                chromGFF.append([chrom,'bin_%s' % (str(ticker)),'',i,i+binSize,'','.','',''])
                ticker+=1
            if not singleGFF:
                unParseTable(chromGFF,gffFolder+'%s_%s_BIN_%s.gff' % (string.upper(genome),string.upper(chrom),str(binSize)),'\t')
            if singleGFF:
                masterGFF+=chromGFF

    if singleGFF:
        unParseTable(masterGFF,gffFolder+'%s_BIN_%s.gff' % (string.upper(genome),str(binSize)),'\t')




    




#==========================================================================
#===================MAKING GFFS OF ENHANCER REGIONS========================
#==========================================================================


def makeEnhancerGFFs(dataFile,gffName,namesList,annotFile,gffFolder,enrichedFolder,window=2000,macs=True):
    '''
    find all possible enhancers.
    enhancers defined as h3k27ac binding sites +/-5kb outside of promoters
    we define center of enhancer as the center of the bound region
    '''

    dataDict = loadDataTable(dataFile)


    if enrichedFolder[-1] != '/':
        enrichedFolder+='/'

    if gffFolder[-1] != '/':
        gffFolder+='/'

    #nameList = ['H128_H3K27AC','H2171_H3K27AC','MM1S_H3K27AC_DMSO','MM1S_H3K27AC_JQ1','U87_H3K27AC','P493-6_T0_H3K27AC','P493-6_T1_H3K27AC','P493-6_T24_H3K27AC']

    #first make the tss collection
    tssGFF = makeTSSGFF(annotFile,5000,5000)
    tssCollection = gffToLocusCollection(tssGFF)

    #make a blank collection to load enhancers into
    enhancerCollection = LocusCollection([],500)

    #don't allow overlapping enhancers
    species = string.upper(dataDict[namesList[0]]['genome'])


    for name in namesList:


        print('finding enhancers in %s' % (name))
        if macs:
            boundCollection = importBoundRegion(enrichedFolder + dataDict[name]['enrichedMacs'],name)
        else:
            boundCollection = importBoundRegion(enrichedFolder + dataDict[name]['enriched'],name)
        for locus in boundCollection.getLoci():
            #make sure an overlapping enhancer doesn't already exist
            if len(tssCollection.getOverlap(locus,'both')) == 0 and len(enhancerCollection.getOverlap(locus,'both')) == 0:
                center = (locus.start()+locus.end())/2
                gffLocus = Locus(locus.chr(),center-window,center+window,'.',locus.ID())
                enhancerCollection.append(gffLocus)


    enhancerGFF = locusCollectionToGFF(enhancerCollection)
    print('Found %s enhancers in %s' % (len(enhancerGFF),gffName))
    unParseTable(enhancerGFF,gffFolder+'%s_ENHANCERS_%s_-%s_+%s.gff' % (species,gffName,window,window),'\t')


#==========================================================================
#===================MAKING GFFS OF ENRICHED REGIONS========================
#==========================================================================


def makeEnrichedGFFs(dataFile,namesList,gffFolder,enrichedFolder,macs=True,window=0):
    '''
    make gffs from enriched regions +/- a window
    '''

    dataDict = loadDataTable(dataFile)
    

    if enrichedFolder[-1] != '/':
        enrichedFolder+='/'

    if gffFolder[-1] != '/':
        gffFolder+='/'


    if len(namesList) == 0:
        namesList = dataDict.keys()
    species = string.upper(dataDict[namesList[0]]['genome'])
    for name in namesList:

        
        print('making enriched region gff in %s' % (name))
        if macs:
            if len(dataDict[name]['enrichedMacs']) ==0 or dataDict[name]['enrichedMacs'] =='NONE':
                continue
            boundCollection = importBoundRegion(enrichedFolder + dataDict[name]['enrichedMacs'],name)
        else:
            if len(dataDict[name]['enriched']) ==0 or dataDict[name]['enriched'] =='NONE':
                continue
            boundCollection = importBoundRegion(enrichedFolder + dataDict[name]['enriched'],name)
        
        boundLoci = boundCollection.getLoci()

        enrichedCollection = LocusCollection([],500)
        for locus in boundLoci:
            searchLocus = makeSearchLocus(locus,int(window),int(window))
            enrichedCollection.append(searchLocus)
            
        enrichedGFF = locusCollectionToGFF(enrichedCollection)
        print('Found %s enriched regions in %s' % (len(enrichedGFF),name))
        unParseTable(enrichedGFF,gffFolder+'%s_ENRICHED_%s_-%s_+%s.gff' % (species,name,window,window),'\t')


#==========================================================================
#===================MAKING GFFS OF PROMOTER REGIONS========================
#==========================================================================

def makePromoterGFF(dataFile,annotFile,promoterFactor,enrichedFolder,gffFolder,window=0,transcribedGeneFile=''):

    '''
    uses a promoter associated factor to define promoter regsion.  Can include a window to extend promoter regions as well as a transcribed gene list to restrict set of genes
    '''
    window = int(window)
    #loading the dataTable
    dataDict = loadDataTable(dataFile)

    #finding the promoter factor in the enriched folder
    formatFolder(enrichedFolder,True)
    formatFolder(gffFolder,True)

    #establishing the output filename
    genome = dataDict[promoterFactor]['genome']
    
    output = '%s%s_PROMOTER_%s_-%s_+%s.gff' % (gffFolder,string.upper(genome),promoterFactor,window,window)

    #getting the promoter factor
    enrichedCollection = importBoundRegion(enrichedFolder + dataDict[promoterFactor]['enrichedMacs'],promoterFactor)

    #making the start dict
    startDict = makeStartDict(annotFile)

    #getting list of transcribed genes
    if len(transcribedGeneFile) > 0:
        transcribedTable = parseTable(transcribedGeneFile,'\t')

        geneList = [line[1] for line in transcribedTable]

    else:
        geneList = startDict.keys()

        
    #now make collection of all transcribed gene TSSs

    tssLoci = []

    for geneID in geneList:
        tssLocus = Locus(startDict[geneID]['chr'],startDict[geneID]['start'][0],startDict[geneID]['start'][0]+1,'.',geneID)
        tssLoci.append(tssLocus)

    tssCollection = LocusCollection(tssLoci,50)

    #a promoter is a single promoter associated enriched region 
    #site that overlaps at most 2 unique genes

    promoterGFF = []
    promoterLoci = enrichedCollection.getLoci()

    for locus in promoterLoci:

        overlappingTSSLoci = tssCollection.getOverlap(locus,'both')
        if len(overlappingTSSLoci) == 0:
            continue
        else:


            geneNames = [startDict[x.ID()]['name'] for x in overlappingTSSLoci]
            geneNames = uniquify(geneNames)

            if len(geneNames) <= 2:
                refseqString = string.join([x.ID() for x in overlappingTSSLoci],',')
                chrom = locus.chr()
                start = locus.start()-window
                end = locus.end()+window
                strand = locus.sense()
                promoterGFF.append([chrom,locus.ID(),'',start,end,'',strand,'',refseqString])

    unParseTable(promoterGFF,output,'\t')



#-------------------------------------------------------------------------#
#                                                                         #
#                             MAPPING TOOLS                               #
#                                                                         #
#-------------------------------------------------------------------------#



#==========================================================================
#===================MAP ENRICHED REGIONS TO GFF============================
#==========================================================================

def mapEnrichedToGFF(dataFile,setName,gffList,cellTypeList,enrichedFolder,mappedFolder,macs=True,namesList=[],useBackground=True):

    '''
    maps enriched regions from a set of cell types to a set of gffs
    tries to make a new folder for each gff
    '''

    dataDict = loadDataTable(dataFile)

    formatFolder(enrichedFolder,True)
    formatFolder(mappedFolder,True)


    for gffFile in gffList:

        gffName = gffFile.split('/')[-1].split('.')[0]
        print('making enriched regions to %s' % (gffName))
        #make the gff into a collection


        try:
            foo = os.listdir(mappedFolder)
        except OSError:
            print('making directory %s to hold mapped enriched files' % (mappedFolder))
            os.system('mkdir %s' % (mappedFolder))


        
        outdir = mappedFolder+gffName+'/'

        try:
            foo = os.listdir(mappedFolder+gffName)
        except OSError:
            os.system('mkdir %s' % (outdir))
        
        #first filter the name list
        cellTypeNameList =[] 
        if len(namesList) == 0:
            namesList = dataDict.keys()
        for name in namesList:
            #check to make sure in the right celltype
            #also make sure to not process WCEs
            if useBackground and dataDict[name]['background'] == 'NONE':
                continue
            cellName = name.split('_')[0]
            if macs == True:
                print(name)
                print(dataDict[name]['enrichedMacs'])
                if cellTypeList.count(cellName) == 1 and dataDict[name]['enrichedMacs'] != 'NONE':
                    cellTypeNameList.append(name)

            else:
                if cellTypeList.count(cellName) == 1 and dataDict[name]['enriched'] != 'NONE':
                    cellTypeNameList.append(name)

        cellTypeNameList.sort()

        mappedGFF = [['GFF_LINE','ID'] + cellTypeNameList]
        #now we go through the gff and fill in stuff
        gffTable = parseTable(gffFile,'\t')

        gffLoci = []
        #making the header line
        for line in gffTable:
            gffLocus = Locus(line[0],line[3],line[4],line[6],line[8])
            gffLine = gffLocus.__str__()
            gffID = line[1]
            
            gffLoci.append(gffLocus)
            mappedGFF.append([gffLine,gffID])
            
        for name in cellTypeNameList:
            print('dataset %s' % (name))
            if macs:
                enrichedCollection = importBoundRegion(enrichedFolder + dataDict[name]['enrichedMacs'],name)
            else:
                enrichedCollection = importBoundRegion(enrichedFolder + dataDict[name]['enriched'],name)
            for i in range(len(gffLoci)):
                if len(enrichedCollection.getOverlap(gffLoci[i],'both')) > 0:
                    mappedGFF[i+1].append(1)
                else:
                    mappedGFF[i+1].append(0)

    out_path = outdir+gffName+'_'+setName+'.txt'
    unParseTable(mappedGFF,out_path,'\t')
    return out_path



#==========================================================================
#===================MAPPING BAMS TO GFFS===================================
#==========================================================================

def mapBams(dataFile,cellTypeList,gffList,mappedFolder,nBin = 200,overWrite =False,rpm=True,nameList = [],extension=200):
    
    '''
    for each gff maps all of the data and writes to a specific folder named after the gff
    can map either by cell type or by a specific name list
    '''

    dataDict = loadDataTable(dataFile)
    if mappedFolder[-1] != '/':
        mappedFolder+='/'
    ticker = 0
    for gffFile in gffList:
        
        #check to make sure gff exists
        try:
            foo = open(gffFile,'r')
        except IOError:
            print('ERROR: GFF FILE %s DOES NOT EXIST' % (gffFile))

        gffName = gffFile.split('/')[-1].split('.')[0]

        #see if the parent directory exists, if not make it
        try:
            foo = os.listdir(mappedFolder)
        except OSError:
            print('%s directory not found for mapped bams. making it' % (mappedFolder))
            os.system('mkdir %s' % (mappedFolder))


        
        #make this directory specifically for the gff
        try:
            foo = os.listdir(mappedFolder+gffName)
        except OSError:
            print('%s directory not found for this gff: %s. making it' % (mappedFolder+gffName,gffName))
            os.system('mkdir %s%s' % (mappedFolder,gffName))

        outdir = mappedFolder+gffName+'/'

        if len(nameList) == 0:
            nameList = dataDict.keys()
        
        

        for name in nameList:
            print ('mapping %s to %s' % (name,gffFile))
            #filter based on celltype
            cellName = name.split('_')[0]
            if cellTypeList.count(cellName) != 1:
                print("this guy didn't get mapped %s" % (name))
                continue
            fullBamFile = dataDict[name]['bam']
            outFile = outdir+gffName+'_'+name+'.gff'

            

            if overWrite:
                cmd1 = "python %s/bamToGFF_turbo.py -e %s -m %s -b %s -i %s -o %s" % (whereAmI,extension,nBin,fullBamFile,gffFile,outFile)
                if rpm:
                    cmd1 += ' -r'
                #cmd1 += ' &'
                print(cmd1)
                os.system(cmd1)

            else:
                try:
                    Foo = open(outFile,'r')
                    print('File %s Already Exists, not mapping' % (outFile))
                except IOError:
                    cmd1 = "python %s/bamToGFF_turbo.py -e %s -m %s -b %s -i %s -o %s" % (whereAmI,extension,nBin,fullBamFile,gffFile,outFile)
                    if rpm:
                        cmd1 += ' -r'
                    #cmd1 += ' &'

                    print(cmd1)
                    os.system(cmd1)




def mapBamsQsub(dataFile,cellTypeList,gffList,mappedFolder,nBin = 200,overWrite =False,nameList = []):
    
    '''
    for each gff maps all of the data and writes to a specific folder named after the gff
    can map either by cell type or by a specific name list
    '''

    dataDict = loadDataTable(dataFile)
    if mappedFolder[-1] != '/':
        mappedFolder+='/'
    #a timestamp to name this pipeline batch of files
    timeStamp =datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    #a random integer ticker to help name files
    randTicker = random.randint(0,10000)
    bashFileName = '/mnt/d0-0/share/bradnerlab/src/cl512/temp/mapBam_%s_%s.sh' % (timeStamp,randTicker)
    bashFile = open(bashFileName,'w')
    ticker = 0
    for gffFile in gffList:
        
        #check to make sure gff exists
        try:
            foo = open(gffFile,'r')
        except IOError:
            print('ERROR: GFF FILE %s DOES NOT EXIST' % (gffFile))

        gffName = gffFile.split('/')[-1].split('.')[0]

        #see if the parent directory exists, if not make it
        try:
            foo = os.listdir(mappedFolder)
        except OSError:
            print('%s directory not found for mapped bams. making it' % (mappedFolder))
            os.system('mkdir %s' % (mappedFolder))


        
        #make this directory specifically for the gff
        try:
            foo = os.listdir(mappedFolder+gffName)
        except OSError:
            print('%s directory not found for this gff: %s. making it' % (mappedFolder+gffName,gffName))
            os.system('mkdir %s%s' % (mappedFolder,gffName))

        outdir = mappedFolder+gffName+'/'

        if len(nameList) == 0:
            nameList = dataDict.keys()
        
        

        for name in nameList:
            print ('mapping %s to %s' % (name,gffFile))
            #filter based on celltype
            cellName = name.split('_')[0]
            if cellTypeList.count(cellName) != 1:
                continue
            fullBamFile = dataDict[name]['bam']
            outFile = outdir+gffName+'_'+name+'.gff'



            if overWrite:
                cmd1 = "python /mnt/d0-0/share/bradnerlab/src/cl512/pipeline/bamToGFF_turbo.py -e 200 -r -m %s -b %s -i %s -o %s &" % (nBin,fullBamFile,gffFile,outFile)
                bashFile.write(cmd1)
                bashFile.write('\n')

            else:
                try:
                    Foo = open(outFile,'r')
                    print('File %s Already Exists, not mapping' % (outFile))
                except IOError:
                    cmd1 = "python /mnt/d0-0/share/bradnerlab/src/cl512/pipeline/bamToGFF_turbo.py -e 200 -r -m %s -b %s -i %s -o %s &" % (nBin,fullBamFile,gffFile,outFile)
                    bashFile.write(cmd1)
                    bashFile.write('\n')



    bashFile.close()
            
    bashCommand = 'qsub %s' % (bashFileName) 
    print(bashCommand)
            #os.system(bashCommand)






def mapBamsBatch(dataFile,gffList,mappedFolder,overWrite =False,namesList = [],extension=200,rpm=True):
    
    '''
    for each gff maps all of the data and writes to a specific folder named after the gff
    can map either by cell type or by a specific name list
    uses bamliquidatorbatch
    '''

    dataDict = loadDataTable(dataFile)
    if mappedFolder[-1] != '/':
        mappedFolder+='/'
    ticker = 0
    for gffFile in gffList:
        
        #check to make sure gff exists
        try:
            foo = open(gffFile,'r')
        except IOError:
            print('ERROR: GFF FILE %s DOES NOT EXIST' % (gffFile))
            sys.exit()

        gffName = gffFile.split('/')[-1].split('.')[0]

        #see if the parent directory exists, if not make it
        mappedFolder = formatFolder(mappedFolder,True)
        outdirRoot = formatFolder(mappedFolder+gffName,True)

      
        if len(namesList) == 0:
            namesList = dataDict.keys()
                

        for name in namesList:
            print ('mapping %s to %s' % (name,gffFile))
            #filter based on celltype
            fullBamFile = dataDict[name]['bam']
        
            #output for the bamliquidator command
            outdir = formatFolder(outdirRoot+name,True)  
            outMatrixFile = outdir+'matrix.txt'
            
            if overWrite:
                mapCmd = bamliquidator_path + " --sense . -e %s --match_bamToGFF -r %s -o %s %s &" % (extension,gffFile, outdir, fullBamFile)    
                print(mapCmd)
                os.system(mapCmd)

            else:
                try:
                    print('checking for outfile %s' % (outMatrixFile))
                    Foo = open(outMatrixFile,'r')
                    print('File %s Already Exists, not mapping' % (outMatrixFile))
                except IOError:
                    mapCmd = bamliquidator_path + " --sense . -e %s --match_bamToGFF -r %s -o %s %s &" % (extension,gffFile, outdir, fullBamFile)                
                    print(mapCmd)
                    os.system(mapCmd)

    #time.sleep(10) #wait 10 seconds before checking for output
    #now initiate another giant loop to check for output and rename it
    for gffFile in gffList:
        
        #check to make sure gff exists
        try:
            foo = open(gffFile,'r')
        except IOError:
            print('ERROR: GFF FILE %s DOES NOT EXIST' % (gffFile))
            sys.exit()

        gffName = gffFile.split('/')[-1].split('.')[0]

        #see if the parent directory exists, if not make it
        mappedFolder = formatFolder(mappedFolder,True)
        #the first outdir of the mapping
        outdirRoot = formatFolder(mappedFolder+gffName,True)

        if len(namesList) == 0:
            namesList = dataDict.keys()
                
        for name in namesList:
            print ('Checking output of %s mapping to %s' % (name,gffFile))
            #filter based on celltype
            fullBamFile = dataDict[name]['bam']
            
            outdir = formatFolder(outdirRoot+name,True)
            matrixFile = outdir + 'matrix.txt'
            countsFile = outdir + 'counts.h5'
            
            #what we want the eventual outfile to look like
            outMatrixFile = outdirRoot+gffName+'_'+name+'.txt'
            outCountsFile = outdirRoot+gffName+'_'+name+'.h5'

            #now make sure the matrix file exists
            try:
                Foo = open(outMatrixFile,'r')    
            except IOError:

                if checkOutput(matrixFile,0.1,2):
                    mvCmd = 'mv %s %s &' % (matrixFile,outMatrixFile)
                    print("Renaming output %s as %s" % (matrixFile,outMatrixFile))
                    os.system(mvCmd)
                else:
                    print("ERROR: No output found for %s mapping to %s" % (name,gffFile))




#==========================================================================
#===================FORMATTING MAPPING SIGNAL==============================
#==========================================================================

def makeSignalTable(dataFile,gffFile,mappedFolder,namesList = [],medianNorm=False,output =''):

    '''
    for each sample, make a dictionary keyed by locus ID
    '''

    dataDict = loadDataTable(dataFile)

    if len(namesList) == 0:
        namesList = dataDict.keys()

    signalDict = {}
    for name in namesList:
        signalDict[name] = defaultdict(float)


    #now start filling in the signal dict
    gffName = gffFile.split('/')[-1].split('.')[0]
    print(gffName)
    for name in namesList:

        print("MAKING SIGNAL DICT FOR %s" % (name))
        
        #try opening the batch mapping output first
        mappedFile = '%s%s/%s_%s.txt' % (mappedFolder,gffName,gffName,name)
        
        if checkOutput(mappedFile,.02,.02):
            print('FOUND MAPPED FILE FOR %s AT %s' % (name,mappedFile))
        else:   
            mappedFile = '%s%s/%s_%s.txt' % (mappedFolder,gffName,gffName,name)
        
        if checkOutput(mappedFile,.02,.02):
            print('FOUND MAPPED FILE FOR %s AT %s' % (name,mappedFile))
        else:
            print('ERROR NO MAPPED FILE FOUND FOR %s' % (name))
            sys.exit()
            
        mappedTable = parseTable(mappedFile,'\t')
        if medianNorm == True:
            medianSignal = numpy.median([float(line[2]) for line in mappedTable[1:]])
        else:
            medianSignal = 1
        
        
        for line in mappedTable[1:]:

            signalDict[name][line[1]] = float(line[2])/medianSignal

    #now make the signal table
    signalTable = []
    header = ['GENE_ID','locusLine'] + namesList
    signalTable.append(header)

    for line in mappedTable[1:]:
        locusID = line[1]
        sigLine = line[0:2] + [signalDict[name][locusID] for name in namesList]
        signalTable.append(sigLine)

    if len(output) == 0:
        return signalTable
    else:
        unParseTable(signalTable,output,'\t')
        return signalTable

#==========================================================================
#============================MAPPING WRAPPER===============================
#==========================================================================

def map_regions(dataFile,gffList,mappedFolder,signalFolder,names_list=[],medianNorm=False,output='',extendReadsTo=200):

    '''
    making a normalized binding signal table at all regions
    '''
    
    #set up a list to return all signal tables made
    signal_table_list = []
    #since each bam has different read lengths, important to carefully normalize quantification
    dataDict = loadDataTable(dataFile)
    dataFile_name = dataFile.split('/')[-1].split('.')[0]

    if len(names_list) == 0:
        names_list = dataDict.keys()
    names_list.sort()
    
    for name in names_list:
        bam = Bam(dataDict[name]['bam'])
        read_length = bam.getReadLengths()[0]
        if int(extendReadsTo) < read_length:
            print('Error: desired overall read extension %s is less than read length %s' % (extendReadsTo,read_length))
            sys.exit()
        bam_extension = int(extendReadsTo)-read_length
        print('For dataset %s using an extension of %s' % (name,bam_extension))
        mapBamsBatch(dataFile,gffList,mappedFolder,overWrite =False,namesList = [name],extension=bam_extension,rpm=True)

    #want a signal table of all datasets to each gff
    print('Writing signal tables for each gff:')
    for gffFile in gffList:
        gffName = gffFile.split('/')[-1].split('.')[0]
        if len(output) == 0:
            signal_table_path = '%s%s_%s_SIGNAL.txt' % (signalFolder,gffName,dataFile_name)
        else:
            signal_table_path = output
        print(signal_table_path)
        makeSignalTable(dataFile,gffFile,mappedFolder,names_list,medianNorm,output =signal_table_path)
        signal_table_list.append(signal_table_path)

    return signal_table_list

#==========================================================================
#===================MAKING GFF LISTS=======================================
#==========================================================================

def makeGFFListFile(mappedEnrichedFile,setList,output,annotFile=''):

    '''
    #setList defines the dataset names to be used
    AND operators within lists, OR operators outside of lists
    [[A,B],[C,D]] = (A AND B) OR (C AND D) for this row
    [[A],[B],[C],[D]] = A OR B OR C OR D for this row

    '''
    if len(annotFile) > 0:
        startDict = makeStartDict(annotFile)

    geneListFile = []
    boundGFFTable = parseTable(mappedEnrichedFile,'\t')
    header = boundGFFTable[0]

    outputFolder = string.join(output.split('/')[0:-1],'/')+'/'
    formatFolder(outputFolder,True)
    
    #convert the setList into column numbers
    columnSet = []
    for bindingSet in setList:
        try:
            columnSet.append([header.index(x) for x in bindingSet])
        except ValueError:
            print('ERROR: not all datasets in binding table')
            exit()

    for i in range(1,len(boundGFFTable),1):
        
        line = boundGFFTable[i]
        refID = line[1]

        #print(columnSet)
        #print(i)
        #if any of these end up being true, the line gets added
        for andColumns in columnSet:
            
            bindingVector = [int(line[x]) for x in andColumns]
            if refID == "NM_133941":
                print(bindingVector)
            #print(bindingVector)
            if bindingVector.count(1) == len(bindingVector):
                if len(annotFile) >0 :
                    geneListFile.append([i,boundGFFTable[i][1],startDict[boundGFFTable[i][1]]['name']])
                else:
                    geneListFile.append([i,boundGFFTable[i][1]])
                break
    print(len(geneListFile))
    unParseTable(geneListFile,output,'\t')

#-------------------------------------------------------------------------#
#                                                                         #
#                            PLOTTING TOOLS                               #
#                                                                         #
#-------------------------------------------------------------------------#

                

#==========================================================================
#===================PLOTTING INDIVIDUAL GENES=============================
#==========================================================================



def callGenePlot(dataFile,geneID,plotName,annotFile,namesList,outputFolder,region='TXN',yScale = 'UNIFORM',plotType='MULTIPLE',nameString =''):

    '''
    calls bamPlot to plot a gene for a given set of bams
    currently works for mm9 and hg18

    '''
    formatFolder(outputFolder,True)

    dataDict = loadDataTable(dataFile)

    startDict= makeStartDict(annotFile)

    if len(geneID) != 0 and startDict.has_key(geneID) == True:
        start = startDict[geneID]['start'][0]
        end = startDict[geneID]['end'][0]
        chrom = startDict[geneID]['chr']
        sense = startDict[geneID]['sense']
        geneLength = abs(end-start)

    if ['TSS','TXN'].count(region) == 1 and len(geneID) == 0:
        print('ERROR: If no gene specified under refseqID, you must enter coordinates in the region field. e.g. chrN:+:1-1000')

    if region == 'TSS':
        
        locusString = '%s:%s:%s-%s' % (chrom,sense,start-5000,start+5000)

    elif region == 'TXN':
        offset =  int(round(.5*geneLength,-3))
        if sense == '+':
            locusString = '%s:%s:%s-%s' % (chrom,sense,start-offset,end+offset)
        else:
            locusString = '%s:%s:%s-%s' % (chrom,sense,end-offset,start+offset)            
    else:
        locusString = region

    if string.upper(plotType) !='MERGE' or len(nameString) == 0:
        nameString = string.join(namesList,',')

    if startDict.has_key(geneID):
        titleString = startDict[geneID]['name'] +'_' + plotName
    else:
        titleString = geneID+'_'+plotName
    colorString = string.join([dataDict[name]['color'] for name in namesList],':')
    bamList = [dataDict[name]['bam'] for name in namesList]
    bamString = string.join(bamList,',')
    genome = string.lower(dataDict[namesList[0]]['genome'])
    os.chdir(whereAmI)
    cmd = "python %s/bamPlot_turbo.py -n %s -t %s -c %s -g %s -p %s -y %s -b %s -i %s -o %s -r --save-temp &" % (whereAmI,nameString,titleString,colorString,genome,plotType,yScale,bamString,locusString,outputFolder)

    #cmd = "python /nfs/young_ata/scripts/bamPlot.py -n %s -t %s -c %s -g hg18 -p multiple -y uniform -b %s -i %s -o %s" % (nameString,titleString,colorString,bamString,locusString,outputFolder)
    print(cmd)
    
    os.system(cmd)

#==========================================================================
#========================BATCH PLOTTING REGIONS============================
#==========================================================================

def callBatchPlot(dataFile,inputFile,plotName,outputFolder,namesList=[],uniform=True,bed ='',plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '',scaleFactorString =''):

    '''
    batch plots all regions in a gff
    if using scale factor, provide the multiplicative scale factors which is 1/rxMMR
    '''
    plotType = string.upper(plotType)
    dataDict = loadDataTable(dataFile)

    if len(namesList) == 0:
        namesList = dataDict.keys()
    
    #we now need to generate all the arguments for bamPlot_turbo.py

    genomeList = [string.lower(dataDict[x]['genome']) for x in namesList]
    if len(uniquify(genomeList)) != 1:
        print("ERROR: CANNOT PLOT DATA FROM MULTIPLE GENOMES")
        sys.exit()
    else:
        genome = genomeList[0]

    #next get the bam string
    bamList = [dataDict[name]['bam'] for name in namesList]
    bamString = string.join(bamList,',')
    
    #inputGFF
    try:
        foo = open(inputFile,'r')
        foo.close()
    except IOError:
        print("ERROR: INPUT FILE NOT READABLE")
        sys.exit()

    #establish the output folder
    outputFolder = formatFolder(outputFolder,True)

    #get the color string
    colorString = string.join([dataDict[name]['color'] for name in namesList],':')

    #get the namesList
    if string.upper(plotType) !='MERGE' or len(nameString) == 0:
        nameString = string.join(namesList,',')
    
    #yScale setting
    if uniform == True:
        yScale = 'UNIFORM'
    else:
        yScale = 'RELATIVE'

    #get the title
    title = plotName

    #figure out if a ChIP-Rx scaling factor is required
    if len(rxGenome) >0:

        #we need to make a list of scaling factors
        scaleFactorList = []
        for name in namesList:
            uniqueID = dataDict[name]['uniqueID']

            #get the rx  bam
            rxBamFile = getTONYInfo(uniqueID,47)
            print('using bamfile %s to scale %s' % (rxBamFile.split('::')[-1],name))
            rxBam = Bam(rxBamFile.split('::')[-1])
            rxMMR = float(rxBam.getTotalReads())/1000000
            scaleFactor = round(1/rxMMR,4)
            print('using a scale factor of %s for %s' % (scaleFactor,name))
            scaleFactorList.append(str(scaleFactor))

        scaleFactorString = string.join(scaleFactorList,',')

    cmd = 'python %sbamPlot_turbo.py -g %s -e %s -b %s -i %s -o %s -c %s -n %s -y %s -t %s -p %s' % (pipelineFolder,genome,extension,bamString,inputFile,outputFolder,colorString,nameString,yScale,title,plotType)
    #scale for RPM
    if rpm == True:
        cmd += ' -r'
    if len(scaleFactorString) > 0:
        if rpm == True:
            print('warning, rpm flag and scale factoring do not mix well')
        cmd += ' --scale %s' % (scaleFactorString)

    if len(bed) > 0:
        cmd += ' --bed %s' % (bed)
    if multiPage:
        cmd += ' --multi-page'
    if debug:
        cmd += ' --save-temp'
    #cmd += ' &'

    print(cmd)
    os.system(cmd)
    return cmd




#=======================================================================
#=======================Plotting TF Motif Corr Heat Maps================
#=======================================================================

def plotCRCCorrMaps(analysis_name,motifBedDir,tf_list_path='',window=50):
    list_files=os.listdir(motifBedDir)
    allLoci=[]
    motif_beds=[]
    for f in list_files:
        fEnd=f.split('.')[-1]
        if fEnd == 'bed':
            motif_beds.append(f)

    print(len(motif_beds))

    temp_dir = '%stmp/' % (motifBedDir)
    figures_dir = '%sfigures/' % (motifBedDir)
    tables_dir = '%stables/' % (motifBedDir)
    
    #making folders
    folderList = [temp_dir,figures_dir,tables_dir]

    for folder in folderList:
        formatFolder(folder,True)

    print(tf_list_path)
    print(len(tf_list_path))

    tf_list=[]
    if len(tf_list_path)>0:
        tf_table=utils.parseTable(tf_list_path,'\t')
        for tf in tf_table:
            print(tf[0])
            tf_list.append(tf[0].upper())
        
        print(tf_list)

        tf_beds=[]
        for bed in motif_beds:
            tf_name=bed.split('_')[0]
            if tf_name in tf_list:
                tf_beds.append(bed)
                print(tf_beds)
            
        motif_beds=tf_beds



    #remove 'track' line from crc bed files and write to tmp file for Rscript use
    for bed in motif_beds:
        bed_path = '%s%s' % (motifBedDir,bed)
        bed_table = utils.parseTable(bed_path,'\t')
        new_table=[]
        for line in bed_table[1:]:
            new_table.append(line)
        tmp_path = '%s%s' % (temp_dir,bed)
        if len(new_table) > 0:
            utils.unParseTable(new_table,tmp_path,'\t')

    for bed in motif_beds:
        TF_name = bed.split('_')[0]
        collection = utils.importBoundRegion('%s%s' %(motifBedDir,bed),TF_name)

        allLoci += collection.getLoci()


    #make stitched bed file
    giant_collection = utils.LocusCollection(allLoci,50)
    stitched_collection = giant_collection.stitchCollection(stitchWindow=50)
    new_bed = utils.locusCollectionToBed(stitched_collection)
    utils.unParseTable(new_bed,'%s%s%s_stitched_bed.bed' % (temp_dir,str(window),analysis_name),'\t')
    utils.unParseTable(new_bed,'%s%s%s_stitched_bed.bed' % (tables_dir,str(window),analysis_name),'\t')


    #Call Rscript to create heatmap figures

    cmd='/storage/cylin/anaconda3/envs/r32_py27/bin/Rscript %splotMotifCorrMatrix.R %s %s %s' % (pipelineFolder,motifBedDir,analysis_name,str(window))
    os.system(cmd)


#-------------------------------------------------------------------------#
#                                                                         #
#                              META TOOLS                                 #
#                                                                         #
#-------------------------------------------------------------------------#


#==========================================================================
#===================MAKING META GFFS OF TXN REGIONS========================
#==========================================================================
        
def makeMetaGFFs(annotFile,gffFolder,genome,geneListFile =''):
    '''
    makes gffs of txn regions for meta genes
    '''
    formatFolder(gffFolder,True)

    genome = string.upper(genome)

        
    startDict = makeStartDict(annotFile)
    
    tssGFF = []
    txnGFF = []
    ttrGFF = []
    
    if len(geneListFile) == 0:
        geneList = startDict.keys()
        geneListName = 'ALL'
    else:
        geneListTable = parseTable(geneListFile,'\t')
        if len(geneListTable[0]) == 1:
            geneList = [line[0] for line in geneListTable]
        else:
            geneList = [line[1] for line in geneListTable]
        geneListName = geneListFile.split('/')[-1].split('.')[0]
    for gene in geneList:

        chrom = startDict[gene]['chr']
        sense = startDict[gene]['sense']
        start = startDict[gene]['start'][0]
        end = startDict[gene]['end'][0]

        if sense == '+':
            tssGFF.append([chrom,gene,'',start-3000,start,'',sense,'',gene])
            txnGFF.append([chrom,gene,'',start,end,'',sense,'',gene])
            ttrGFF.append([chrom,gene,'',end,end+3000,'',sense,'',gene])
        else:
            tssGFF.append([chrom,gene,'',start,start+3000,'',sense,'',gene])
            txnGFF.append([chrom,gene,'',end,start,'',sense,'',gene])
            ttrGFF.append([chrom,gene,'',end-3000,end,'',sense,'',gene])


    unParseTable(tssGFF,'%s%s_TSS_%s_-3000_+0.gff' % (gffFolder,genome,geneListName),'\t')
    unParseTable(txnGFF,'%s%s_TXN_%s_-0_+0.gff' % (gffFolder,genome,geneListName),'\t')
    unParseTable(ttrGFF,'%s%s_TTR_%s_-0_+3000.gff' % (gffFolder,genome,geneListName),'\t')




#==========================================================================
#=======================MAPPING BAMS FOR METAS=============================
#==========================================================================

def mapMetaBams(dataFile,metaName,gffList,cellTypeList,metaFolder,nameList= [],overwrite=False):

    '''
    calls bam mapping for all pol2 and respective wce samples
    '''

    formatFolder(metaFolder,True)

    dataDict = loadDataTable(dataFile)
    if len(nameList) == 0:
        nameList = dataDict.keys()
    outdir = metaFolder + metaName + '/'
    print('out directory is %s' % outdir)

    try:
        foo = os.listdir(outdir)
    except OSError:
        os.system('mkdir %s' % (outdir))
    
    gffString = string.join(gffList,',')

    timeStamp =datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    #a random integer ticker to help name files
    randTicker = random.randint(0,10000)
    bashFileName = '%smapMetaCall_%s_%s.sh' % (metaFolder,timeStamp,randTicker)
    bashFile = open(bashFileName,'w')


    for cellType in cellTypeList:
        

        cellTypeNames =[]
        
        #find all of the pol2 stuff
        for name in nameList:
            if name.count(cellType) == 1:
                cellTypeNames+=[name]
        print(cellTypeNames)
        for name in cellTypeNames:


                    
            bamFile = dataDict[name]['bam']

            cmd = 'python %s/makeBamMeta.py -c -n %s -g %s -b %s -o %s &' % (whereAmI,name,gffString,bamFile,outdir)
            if overwrite == False:
                try:
                    foo = open('%s%s_metaSettings.txt' % (outdir,name),'r')
                    print('meta for %s already exists, skipping.' % (name))
                    continue
                except IOError:
                    #print(cmd)
                    bashFile.write(cmd)
                    bashFile.write('\n')
                    #subprocess.call(cmd,shell=True)
                    #os.system(cmd)
            else:

                #print(cmd)
                bashFile.write(cmd)
                bashFile.write('\n')

                #subprocess.call(cmd,shell=True)
                #os.system(cmd)
        bashFile.close()
        print(bashFileName)


#==========================================================================
#=========================FINISHING METAS==================================
#==========================================================================

def finishMetas(metaFolder,settingsFileList=[]):

    '''
    manually finishes off meta code in case bad things happened

    '''

    if metaFolder[-1] != '/':
        metaFolder+='/'

    metaFiles = os.listdir(metaFolder)
    if len(settingsFileList) == 0:
        settingsFileList = filter(lambda x:x.count('metaSettings') == 1,metaFiles)

    for settingsFile in settingsFileList:
        settingsName = settingsFile[0:-17]
        cmd = ' python /mnt/d0-0/share/bradnerlab/src/cl512/pipeline/makeBamMeta.py -c -f -n %s -o %s' % (settingsName,metaFolder)
        
        print(cmd)
        os.system(cmd)
    

#==========================================================================
#===================MAKING ORDERED HEATMAPS================================
#==========================================================================

    
def callHeatPlotOrdered(dataFile,gffFile,namesList,orderByName,geneListFile,outputFolder,mappedFolder,relative=False,useBackground=False):

    '''
    calls a heatmap ordered by a single dataset
    will spit out a series of heatmaps all ordered by the same reference dataset
    bound table is the binding table for enriched regions in this gff
    boundList is the OR set of requirements necessary to qualify a 
    '''

    dataDict = loadDataTable(dataFile)
    
    #if a blank geneListFile is given, set it to 'NONE'
    if geneListFile == '':
        geneListFile = 'NONE'

    if mappedFolder[-1] != '/':
        mappedFolder+='/'

    if outputFolder[-1] != '/':
        outputFolder+='/'
        
    formatFolder(outputFolder,True)

    gffName = gffFile.split('/')[-1].split('.')[0]
    #get all of the mappedGFFs

    referenceMappedGFF = mappedFolder + gffName + '/' + gffName + '_'+orderByName + '.gff'
    
    for name in namesList:

        mappedGFF = mappedFolder + gffName + '/' + gffName + '_'+name + '.gff'

        if useBackground:
            backgroundName = dataDict[name]['background']
            if backgroundName == 'NONE':
                backgroundGFF = 'NONE'
            backgroundGFF = mappedFolder + gffName + '/' + gffName + '_' + backgroundName + '.gff'
        else:
            backgroundGFF = 'NONE'
        
        color = dataDict[name]['color']
        output = outputFolder + '%s_%s_%s_order.png' % (gffName,name,orderByName)

        cmd = "Rscript %sheatMapOrdered.R %s %s %s %s %s" % (pipelineFolder,referenceMappedGFF,mappedGFF,color,output,geneListFile)
        if relative:
            cmd += ' 1'
        else:
            cmd += ' 0'

        #now add the background stuff
        cmd += ' %s' % (backgroundGFF)




        print(cmd)
        os.system(cmd)


#==========================================================================
#==============================CALLING ROSE================================
#==========================================================================



def callRose(dataFile,macsEnrichedFolder,parentFolder,namesList=[],extraMap = [],inputFile='',tss=2500,stitch=12500,bashFileName ='',mask=''):

    '''
    calls rose w/ standard parameters
    '''

    dataDict = loadDataTable(dataFile)
    if len(namesList) == 0:
        namesList = dataDict.keys()
    print(namesList)
    #a timestamp to name this pipeline batch of files
    timeStamp =datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    #a random integer ticker to help name files
    randTicker = random.randint(0,10000)

    formatFolder(parentFolder,True)

    if len(bashFileName) == 0:
        bashFileName = '%srose_%s_%s.sh' % (parentFolder,timeStamp,randTicker)
    bashFile = open(bashFileName,'w')
    bashFile.write("cd /ark/home/cl512/src/rose/")
    bashFile.write('\n')

    mapString = [dataDict[name]['bam'] for name in extraMap]
    mapString = string.join(mapString,',')

    for name in namesList:
        #print name
        genome = dataDict[name]['genome']
        bamFile = dataDict[name]['bam']
        
        backgroundName = dataDict[name]['background']
        if dataDict.has_key(backgroundName):
            backgroundBamFile = dataDict[backgroundName]['bam']
            hasBackground = True
        else:
            hasBackground = False

        if len(inputFile) == 0:
            macsFile = "%s%s" % (macsEnrichedFolder,dataDict[name]['enrichedMacs'])
        else:
            macsFile = inputFile
        outputFolder = "%s%s_ROSE" % (parentFolder,name)

        roseCmd = 'python ROSE_main_turbo.py -g %s -i %s -r %s -o %s -t %s -s %s' % (genome,macsFile,bamFile,outputFolder,tss,stitch)

        if hasBackground:
            roseCmd +=' -c %s' % (backgroundBamFile)
        if len(mapString) > 0:
            roseCmd +=' -b %s' % (mapString)
        if len(mask) >0:
            roseCmd += ' --mask %s' % (mask)

        roseCmd += ' &'
        bashFile.write(roseCmd)
        bashFile.write('\n')


    bashFile.close()

    print ('Wrote rose commands to %s' % (bashFileName))
    return bashFileName




def callRose2(dataFile,macsEnrichedFolder,parentFolder,namesList=[],extraMap = [],inputFile='',tss=2500,stitch='',bashFileName ='',mask='',useBackground=True,py27_path =''):

    '''
    calls rose w/ standard parameters
    '''
    if py27_path =='':
        py27_path = 'python'

    dataDict = loadDataTable(dataFile)
    if len(namesList) == 0:
        namesList = dataDict.keys()

    #a timestamp to name this pipeline batch of files
    timeStamp =datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    #a random integer ticker to help name files
    randTicker = random.randint(0,10000)

    formatFolder(parentFolder,True)

    if len(bashFileName) == 0:
        bashFileName = '%srose_%s_%s.sh' % (parentFolder,timeStamp,randTicker)
    bashFile = open(bashFileName,'w')
    bashFile.write("cd %s" % (pipelineFolder))
    bashFile.write('\n')

    mapString = [dataDict[name]['bam'] for name in extraMap]
    mapString = string.join(mapString,',')

    for name in namesList:
        #print name
        genome = dataDict[name]['genome']
        bamFile = dataDict[name]['bam']
        
        backgroundName = dataDict[name]['background']
        if useBackground and dataDict.has_key(backgroundName):
            backgroundBamFile = dataDict[backgroundName]['bam']
            hasBackground = True
        else:
            hasBackground = False

        if len(inputFile) == 0:
            macsFile = "%s%s" % (macsEnrichedFolder,dataDict[name]['enrichedMacs'])
        else:
            macsFile = inputFile
        outputFolder = "%s%s_ROSE" % (parentFolder,name)
        print(name)
        bashFile.write('#running ROSE2 on %s\n' % (name))
        roseCmd = '%s ROSE2_main.py -g %s -i %s -r %s -o %s -t %s' % (py27_path,genome,macsFile,bamFile,outputFolder,tss)

        if len(str(stitch)) > 0:
            roseCmd += ' -s %s' % (stitch)
        if hasBackground:
            roseCmd +=' -c %s' % (backgroundBamFile)
        if len(mapString) > 0:
            roseCmd +=' -b %s' % (mapString)
        if len(mask) >0:
            roseCmd += ' --mask %s' % (mask)

        bashFile.write(roseCmd)
        bashFile.write('\n\n')


    bashFile.close()

    print ('Wrote rose commands to %s' % (bashFileName))
    return bashFileName



def callRose2Slurm(dataFile,macsEnrichedFolder,parentFolder,namesList=[],extraMap = [],inputFile='',tss=2500,stitch='',bashFileName ='',mask='',useBackground=True,projectName=''):

    '''
    calls rose w/ standard parameters
    '''

    #load the data dict
    dataDict = loadDataTable(dataFile)
    if len(namesList) == 0:
        namesList = dataDict.keys()


    #for recording purposes
    #first set the project name
    if len(projectName) == 0:
        #draw from the data file
        projectName = '%s_ROSE2' % (dataFile.split('/')[-1].split('.txt')[0])

    #a timestamp to name this pipeline batch of files
    timeStamp =datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    #a random integer ticker to help name files
    randTicker = random.randint(0,10000)

    formatFolder(parentFolder,True)

    if len(bashFileName) == 0:
        bashFileName = '%srose_%s_%s.sh' % (parentFolder,timeStamp,randTicker)
    bashFile = open(bashFileName,'w')

    ts = time.time()
    timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%Hh%Mm%Ss')
    cmd = '#!/usr/bin/bash'
    bashFile.write(cmd+'\n')
    cmd = '#SBATCH --output=/storage/cylin/grail/slurm_out/ROSE2_%s_%s' % (projectName,timestamp) + '_%j.out # Standard output and error log'
    bashFile.write(cmd+'\n')
    cmd = '#SBATCH -e /storage/cylin/grail/slurm_out/ROSE2_%s_%s' % (projectName,timestamp) + '_%j.err # Standard output and error log'
    bashFile.write(cmd+'\n')
    cmd = '#SBATCH -n 8'
    bashFile.write(cmd+'\n')
    cmd = '#SBATCH --mem 32768'
    bashFile.write(cmd+'\n')


    cmd = 'pwd; hostname; date'
    bashFile.write(cmd+'\n\n\n\n')


    bashFile.write("cd %s" % (whereAmI))
    bashFile.write('\n')

    mapString = [dataDict[name]['bam'] for name in extraMap]
    mapString = string.join(mapString,',')

    for name in namesList:
        #print name
        genome = dataDict[name]['genome']
        bamFile = dataDict[name]['bam']

        backgroundName = dataDict[name]['background']
        if useBackground and dataDict.has_key(backgroundName):
            backgroundBamFile = dataDict[backgroundName]['bam']
            hasBackground = True
        else:
            hasBackground = False

        if len(inputFile) == 0:
            macsFile = "%s%s" % (macsEnrichedFolder,dataDict[name]['enrichedMacs'])
        else:
            macsFile = inputFile
        outputFolder = "%s%s_ROSE" % (parentFolder,name)

        roseCmd = 'python ROSE2_main.py -g %s -i %s -r %s -o %s -t %s' % (genome,macsFile,bamFile,outputFolder,tss)

        if len(str(stitch)) > 0:
            roseCmd += ' -s %s' % (stitch)
        if hasBackground:
            roseCmd +=' -c %s' % (backgroundBamFile)
        if len(mapString) > 0:
            roseCmd +=' -b %s' % (mapString)
        if len(mask) >0:
            roseCmd += ' --mask %s' % (mask)

        roseCmd += ''


        bashFile.write(roseCmd)
        bashFile.write('\n')


    bashFile.close()

    print ('Wrote rose commands to %s' % (bashFileName))
    return bashFileName


#-------------------------------------------------------------------------#
#                                                                         #
#                               CRC TOOLS                                 #
#                                                                         #
#-------------------------------------------------------------------------#



def call_crc(analysis_name,enhancer_path,subpeak_path,activity_path,genome,crc_folder='',args ='',py27_path =''):

    '''
    runs crc
    '''
    
    #if no python path is provide use default python
    if py27_path == '':
        py27_path = 'python'

    if len(crc_folder) == 0:
        
        crc_folder = utils.formatFolder('./crc',True)
    else:
        crc_folder = utils.formatFolder(crc_folder,True)


    output_folder = utils.formatFolder('%s%s' % (crc_folder,analysis_name),True)

    crc_cmd = '%s %scrc/CRC3.py -e %s -g %s -o %s -n %s -s %s -a %s' % (py27_path,pipelineFolder,enhancer_path,genome,output_folder,analysis_name,subpeak_path,activity_path)
    if len(args) > 0:
        if args[0] != ' ': #adds an extra space
            args = ' ' + args

        crc_cmd += args
    crc_bash_path = '%s%s_crc.sh' % (crc_folder,analysis_name)

    crc_bash = open(crc_bash_path,'w')
    crc_bash.write('#!/usr/bin/bash\n\n')

    crc_bash.write('#running crc for %s\n' % (analysis_name))
    crc_bash.write(crc_cmd +'\n')
    crc_bash.close()

    print('wrote crc command for %s to %s' % (analysis_name,crc_bash_path))
    return crc_bash_path





#-------------------------------------------------------------------------#
#                                                                         #
#                          EXPRESSION TOOLS                               #
#                                                                         #
#-------------------------------------------------------------------------#

def mapHisat(dataFile,namesList=[],projectFolder='',useSRA=False,pCount=16,Launch=True):
    
   # '''
   # maps using hisat2 if useSRA is flagged will try to extract an SRA ID from the fastq path and call directly
   # '''
   #load the datasets
    dataDict = loadDataTable(dataFile)

    #loading the datasets we want to process
    if len(namesList) == 0:
        namesList = dataDict.keys()

    #want to write to the bam folder which we can deduce from the genomes of the datasets
    #we want to make sure everyone has the same genome or else this is scary

    genomeList = [dataDict[name]['genome'] for name in dataDict.keys()]

    if len(utils.uniquify(genomeList)) > 1:
        print('OH HECK NO YOU CANT ALIGN MULTIPLE GENOMES THAT WOULD BE STUPID')
        sys.exit()

    genome = string.upper(genomeList[0])

    bamFolder = utils.formatFolder(dataDict[namesList[0]]['folder'])
    hisatIndexDictionary = {'HG38':'/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Hisat2Index/hg38',
                           'MM9': '/storage/cylin/grail/genomes/Mus_musculus/UCSC/mm9/Sequence/Hisat2Index/mm9',
                           'RN6_ERCC': '/storage/cylin/grail/genomes/Rattus_norvegicus/UCSC/rn6/Sequence/Hisat2Index_ERCC/rn6_ercc',
                           'RN6': '/storage/cylin/grail/genomes/Rattus_norvegicus/UCSC/rn6/Sequence/Hisat2Index/rn6',
                           'HG19_ERCC': '/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Hisat2Index_ERCC/hg19_ercc',
                           'HG19': '/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Hisat2Index/hg19', 
                           'MM10': '/storage/cylin/grail/genomes/Mus_musculus/UCSC/mm10/Sequence/Hisat2Index/mm10' 
                          }
    
    hisatIndex = hisatIndexDictionary[genome]


  #now we can loop through the datasets
    for name in namesList:
        bashFilePath = '%s%s_HISAT2' % (projectFolder,name)
        bashFile = open(bashFilePath,'w')
        bashFile.write('#!/usr/bin/bash\n') #shebang line plus end of line characters
        bashFile.write('\n\n\n')
       
        #first write a line to cd into the bam folder
        bashFile.write('cd %s\n' %(bamFolder))
        
        uniqueID = dataDict[name]['uniqueID']
        ts = time.time()
        timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%Hh%Mm%Ss')
        cmd = '#SBATCH --output=/storage/cylin/grail/slurm_out/hisat2_%s_%s' % (name,timestamp) + '_%j.out # Standard output and error log'
        bashFile.write(cmd+'\n')
        cmd = '#SBATCH -e /storage/cylin/grail/slurm_out/hisat2_%s_%s' % (name,timestamp) + '_%j.err # Standard output and error log'
        bashFile.write(cmd+'\n')

        cmd = 'pwd; hostname; date'
        bashFile.write(cmd+'\n')
        bashFile.write('\n\n\n')

        bashFile.write('#===================\n')
        bashFile.write('#PROCESSING %s\n' %(name))
        bashFile.write('echo "processing %s"\n' % (name))

        outputSam = '%s%s.%s.sam' % (bamFolder,uniqueID,genome)
        outputBam = '%s%s.%s.bam' % (bamFolder,uniqueID,genome)
        outputSortedBam = '%s%s.%s.sorted' % (bamFolder,uniqueID,genome)

        if useSRA:
            srrID = dataDict[name]['fastq'].split('/')[-2]
            alignCmd = 'hisat2 -p %s --no-unal -x %s --sra-acc %s -S %s' % (pCount,hisatIndex,srrID,outputSam)
        else:
            #check for paired end
            fastqPath = dataDict[name]['fastq']
            if fastqPath.count('::') == 1:
                #this is paried end
                [fastqPath_1,fastqPath_2] = fastqPath.split('::')
                alignCmd = 'hisat2 -p %s --no-unal -x %s -1 %s -2 %s -S %s' % (pCount,hisatIndex,fastqPath_1,fastqPath_2,outputSam)
            else:
                alignCmd = 'hisat2 -p %s --no-unal -x %s -U %s -S %s' % (pCount,hisatIndex,fastqPath,outputSam)

        bashFile.write(alignCmd)
        bashFile.write('\n')

      # now convert the sam to a bam
        generateBamCmd = '/usr/bin/samtools view -bS %s > %s' % (outputSam,outputBam)
        bashFile.write(generateBamCmd)
        bashFile.write('\n')

      # now we need to sort the bam
        sortBamCmd = '/usr/bin/samtools sort %s %s' % (outputBam,outputSortedBam)
        bashFile.write(sortBamCmd)
        bashFile.write('\n')

      # now we need to index the bam
        indexBamCmd = '/usr/bin/samtools index %s.bam' % (outputSortedBam)
        bashFile.write(indexBamCmd)
        bashFile.write('\n')

      #now we need to delete the sam
        deleteSamCmd = 'rm %s' % (outputSam)
        bashFile.write(deleteSamCmd)
        bashFile.write('\n')

        #now we need to delete the unsorted bam
        deleteBamCmd = 'rm %s' % (outputBam)
        bashFile.write(deleteBamCmd)
        bashFile.write('\n')
        bashFile.write('\n\n\n')
        bashFile.close()
        if Launch:
            time.sleep(1)
            cmd = "sbatch -n %s --mem 32768 %s &" % (pCount,bashFilePath)
            os.system(cmd)




def makeCuffTable(dataFile,analysisName,gtfFile,cufflinksFolder,groupList=[],bashFileName = '',useERCC = True):

    '''
    call cuffquant on each bam individually
    and then string the cbx files into cuffnorm
    groupList = [['A_1','A_2'],['B_1','B_2']]
    '''

    def long_substr(data):
        '''
        helper function to find longest substring for group naming
        '''
        substr = ''
        if len(data) > 1 and len(data[0]) > 0:
            for i in range(len(data[0])):
                for j in range(len(data[0])-i+1):
                    if j > len(substr) and all(data[0][i:i+j] in x for x in data):
                        substr = data[0][i:i+j]
        return substr
    
    dataDict = loadDataTable(dataFile)

    #if no grouplist is given
    #run every dataset as a single group
    #for now assumes that every dataset given is RNA Seq
    if len(groupList) == 0:
        namesList = dataDict.keys()
        namesList.sort()
        groupList = [[x] for x in namesList]
        namesString = ','.join(namesList)

    else:
        #only a single name per group
        namesList =[]
        namesStringList = []
        groupTicker = 1
        for group in groupList:
            
            namesList+=group
            coreName = long_substr(group)
            if len(coreName) ==0:
                coreName = '%s_GROUP_%s' % (analysisName,groupTicker)
            else:
                if '-_.'.count(coreName[-1]) == 1:  #get rid of any separators for a core name
                    coreName = coreName[:-1]
            namesStringList.append(coreName)
        print('Using the following as group names')
        print(namesStringList)
        if len(utils.uniquify(namesStringList)) != len(groupList):
            print('Error: only found %s unique group strings to go with %s groups' % (len(utils.uniquify(namesStringList)),len(groupList)))

            sys.exit()
                      
            groupTicker+=1
        namesString = ','.join(namesStringList)
            
    cufflinksFolder = formatFolder(cufflinksFolder,True)

    #let's do this in bashfile format
    if len(bashFileName) ==0:
        bashFileName = '%scuffquant.sh' % (cufflinksFolder)
        
    
    bashFile = open(bashFileName,'w')

    bashFile.write('#!/usr/bin/bash\n')

    bashFile.write('cd %s\n\n' % (cufflinksFolder))

    bashFile.write("echo 'making cuffquant folders'\n")

    for name in namesList:
        bashFile.write('mkdir %s\n' % (name))

    bashFile.write("\necho 'calling cuffquant'\n")

    cuffquantList = [] # create a list to store cuffquant .cxb outputs so we can check for completeness
    for name in namesList:
        bamFileName = dataDict[name]['bam']
        bashFile.write('cuffquant -p 4 -o %s%s/ %s %s --library-type fr-firststrand\n' % (cufflinksFolder,name,gtfFile,bamFileName))
        cuffquantList.append('%s%s/abundances.cxb' % (cufflinksFolder,name))


    #if we want to have python run this as opposed to making a bash file
    # #check for output
    # for cuffquantFile in cuffquantList:

    #     if checkOutput(cuffquantFile,5,60):
    #         print "FOUND CUFFQUANT OUTPUT FOR %s" % (cuffquantFile)
            
    #     else:
            
        
    #now we want to string together all of the abundances.cxb files to run cuffnorm
    #cuff norm gives you the opportunity to string together replicates
    #gotta figure out the right way to designate sample groups

    cxbList = []
    for group in groupList:
        
        groupString = ','.join(['%s%s/abundances.cxb' % (cufflinksFolder,name) for name in group])
        cxbList.append(groupString)

    cxbString = ' '.join(cxbList)

    #set up the analysis output folders
    cuffnormFolder = formatFolder('%s%s_cuffnorm' % (cufflinksFolder,analysisName),True)
    rOutputFolder = formatFolder('%s%s_cuffnorm/output/' % (cufflinksFolder,analysisName),True)

    #now run the cuffnorm    
    bashFile.write("\necho 'running cuffnorm command'\n")

    
    cuffNormCmd = 'cuffnorm -p 4 -o %s%s_cuffnorm/ -L %s %s %s --library-type fr-firststrand\n' % (cufflinksFolder,analysisName,namesString,gtfFile,cxbString)

    bashFile.write(cuffNormCmd + '\n')


    #now we'll want to pipe the output into the R script for RNA_Seq normalization
    geneFPKMFile = '%s%s_cuffnorm/genes.fpkm_table' % (cufflinksFolder,analysisName)


    if useERCC:
        rCmd = '#Rscript %snormalizeRNASeq.R %s %s %s %s TRUE\n' % (pipelineFolder,geneFPKMFile,rOutputFolder,analysisName,namesString)
    else:
        rCmd = '#Rscript %snormalizeRNASeq.R %s %s %s %s FALSE\n' % (pipelineFolder,geneFPKMFile,rOutputFolder,analysisName,namesString)
    bashFile.write(rCmd)
    bashFile.close()




def makeCuffTableSlurm(dataFile,analysisName,gtfFile,cufflinksFolder,groupList=[],bashFileName = ''):

    '''
    call cuffquant on each bam individually
    and then string the cbx files into cuffnorm
    groupList = [['A_1','A_2'],['B_1','B_2']]
    '''

    def long_substr(data):
        '''
        helper function to find longest substring for group naming
        '''
        substr = ''
        if len(data) > 1 and len(data[0]) > 0:
            for i in range(len(data[0])):
                for j in range(len(data[0])-i+1):
                    if j > len(substr) and all(data[0][i:i+j] in x for x in data):
                        substr = data[0][i:i+j]
        return substr

    dataDict = loadDataTable(dataFile)

    #if no grouplist is given
    #run every dataset as a single group
    #for now assumes that every dataset given is RNA Seq
    if len(groupList) == 0:
        namesList = dataDict.keys()
        namesList.sort()
        groupList = [[x] for x in namesList]
        namesString = ','.join(namesList)

    else:
        #only a single name per group
        namesList =[]
        namesStringList = []
        groupTicker = 1
        for group in groupList:

            namesList+=group
            coreName = long_substr(group)
            if len(coreName) ==0:
                coreName = '%s_GROUP_%s' % (analysisName,groupTicker)
            else:
                if '-_.'.count(coreName[-1]) == 1:  #get rid of any separators for a core name
                    coreName = coreName[:-1]
            namesStringList.append(coreName)
        print('Using the following as group names')
        print(namesStringList)
        if len(utils.uniquify(namesStringList)) != len(groupList):
            print('Error: only found %s unique group strings to go with %s groups' % (len(utils.uniquify(namesStringList)),len(groupList)))
            sys.exit()

            groupTicker+=1
        namesString = ','.join(namesStringList)

    cufflinksFolder = formatFolder(cufflinksFolder,True)

    #let's do this in bashfile format
    if len(bashFileName) ==0:
        bashFileName = '%scuffquant.sh' % (cufflinksFolder)


    bashFile = open(bashFileName,'w')

    bashFile.write('#!/usr/bin/bash\n')

    ts = time.time()
    timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%Hh%Mm%Ss')
    cmd = '#SBATCH --output=/storage/cylin/grail/slurm_out/cufflinks_%s_%s' % (analysisName,timestamp) + '_%j.out # Standard output and error log'
    bashFile.write(cmd+'\n')
    cmd = '#SBATCH -e /storage/cylin/grail/slurm_out/cufflinks_%s_%s' % (analysisName,timestamp) + '_%j.err # Standard output and error log'
    bashFile.write(cmd+'\n')

    cmd = 'pwd; hostname; date'
    bashFile.write(cmd+'\n\n\n\n')


    bashFile.write('cd %s\n\n' % (cufflinksFolder))

    bashFile.write("echo 'making cuffquant folders'\n")

    for name in namesList:
        bashFile.write('mkdir %s\n' % (name))

        bashFile.write("\necho 'calling cuffquant'\n")

    cuffquantList = [] # create a list to store cuffquant .cxb outputs so we can check for completeness
    for name in namesList:
        bamFileName = dataDict[name]['bam']
        bashFile.write('cuffquant -p 4 -o %s%s/ %s %s &\n' % (cufflinksFolder,name,gtfFile,bamFileName))
        cuffquantList.append('%s%s/abundances.cxb' % (cufflinksFolder,name))


    #if we want to have python run this as opposed to making a bash file
    # #check for output
    # for cuffquantFile in cuffquantList:

    #     if checkOutput(cuffquantFile,5,60):
    #         print "FOUND CUFFQUANT OUTPUT FOR %s" % (cuffquantFile)

    #     else:


    #now we want to string together all of the abundances.cxb files to run cuffnorm
    #cuff norm gives you the opportunity to string together replicates
    #gotta figure out the right way to designate sample groups

    cxbList = []
    for group in groupList:

        groupString = ','.join(['%s%s/abundances.cxb' % (cufflinksFolder,name) for name in group])
        cxbList.append(groupString)

    cxbString = ' '.join(cxbList)

    #set up the analysis output folders
    cuffnormFolder = formatFolder('%s%s_cuffnorm' % (cufflinksFolder,analysisName),True)
    rOutputFolder = formatFolder('%s%s_cuffnorm/output/' % (cufflinksFolder,analysisName),True)

    #now run the cuffnorm
    bashFile.write("\necho 'running cuffnorm command'\n")


    cuffNormCmd = 'cuffnorm -p 4 -o %s%s_cuffnorm/ -L %s %s %s\n' % (cufflinksFolder,analysisName,namesString,gtfFile,cxbString)

    bashFile.write(cuffNormCmd + '\n')


    #now we'll want to pipe the output into the R script for RNA_Seq normalization
    geneFPKMFile = '%s%s_cuffnorm/genes.fpkm_table' % (cufflinksFolder,analysisName)




    rCmd = '#Rscript %snormalizeRNASeq.R %s %s %s %s TRUE\n' % (pipelineFolder,geneFPKMFile,rOutputFolder,analysisName,namesString)

    bashFile.write(rCmd)
    bashFile.close()




#-------------------------------------------------------------------------#
#                                                                         #
#                              GECKO TOOLS                                #
#                                                                         #
#-------------------------------------------------------------------------#





def processGecko(dataFile,geckoFolder,namesList = [],overwrite=False,scoringMethod='WtSum'):

    '''
    processes the GECKO bams from the dataFile/namesList
    '''

    dataDict = loadDataTable(dataFile)

    #making the gecko folder
    geckoFolder = formatFolder(geckoFolder,True)

    #if no names are given, process all the data
    if len(namesList) == 0:

        namesList = dataDict.keys()


    for name in namesList:
        
        #skip if a background set
        if string.upper(dataDict[name]['background']) == 'NONE' or string.upper(dataDict[name]['background']) == '':
            continue

        #check for the bam
        try:
            bam = open(dataDict[name]['bam'],'r')
            hasBam = True

        
        except IOError:
            hasBam = False

        #check for background
        try:
            backgroundName = dataDict[name]['background']
            backbroundBam = open(dataDict[backgroundName]['bam'],'r')
            hasBackground = True
        except (IOError, KeyError) as e:
            hasBackground = False

        if not hasBam:
            print('no bam found for %s. GECKO processing not called' % (name))
            continue

        if not hasBackground:
            print('no background bam %s found for dataset %s. GECKO processing not called' % (backgroundName,name))
            continue

        #getting bam names
        testBamFileName = dataDict[name]['bam']
        controlBamFileName = dataDict[backgroundName]['bam']

        #names string
        namesString = '%s,%s' % (name,backgroundName)

        genome = string.lower(dataDict[name]['genome'])
      

        #now set up the process Gecko command

        processGeckoCmd = 'python %sprocessGeckoBam.py -t %s -c %s -g %s -n %s -s %s -o %s' % (pipelineFolder,testBamFileName,controlBamFileName,genome,namesString,scoringMethod,geckoFolder)
        print("RUNNING GECKO PROCESSING FOR %s WITH BACKGROUND %s IN BUILD %s WITH %s SCORING" % (name,backgroundName,genome,scoringMethod))
        print(processGeckoCmd)



        os.system(processGeckoCmd)








#=========================================================================================================
#EVERYTHING BELOW IS EITHER DEPRECATED OR HAS NOT BEEN INTEGRATED INTO PIPELINE YET


#==========================================================================
#=======================PLOTTING GROUPS OF GENES===========================
#==========================================================================

# def plotGeneList(dataFile,annotFile,geneList,namesList,outputFolder,upsearch,downsearch,yScale = 'UNIFORM',byName = True):

#     '''
#     plots the txn region for a bunch of genes in a given window. uses either name or ID
#     deprecating this puppy
#     '''
    
#     startDict = makeStartDict(annotFile)
#     if byName:
#         geneList = [string.upper(x) for x in geneList]
#         nameDict = defaultdict(list)
        
#         for key in startDict.keys():
#             nameDict[startDict[key]['name']].append(key)
            
#         refIDList = []

#         for geneName in geneList:
#             mouseName = string.lower(geneName)
#             mouseName = string.upper(mouseName[0]) + mouseName[1:]
#             #tries both the human and mouse nomenclature
#             refIDs = nameDict[string.upper(geneName)] + nameDict[mouseName]
#             if len(refIDs) ==0:
#                 print('Gene name %s not in annotation file %s' % (geneName,annotFile))
#                 continue
#             refNumbers = [x.split('_')[-1] for x in refIDs]
#             #take the lowest refseq number for the gene name as this is usually the best annotation
#             refIndex = refNumbers.index(min(refNumbers))
#             refIDList.append(refIDs[refIndex])
#             print('Gene %s corresponds to ID %s' % (geneName,refIDs[refIndex]))

#     else:
#         refIDList = geneList

#     for refID in refIDList:

#         chrom = startDict[refID]['chr']
#         sense = startDict[refID]['sense']
#         if sense == '+':
#             start = startDict[refID]['start'][0] - upsearch
#             stop = startDict[refID]['end'][0] + downsearch
#         else:
#             start = startDict[refID]['start'][0] + upsearch
#             stop = startDict[refID]['end'][0] - downsearch

#         geneName = startDict[refID]['name']
#         plotName = '_-%s_+%s' % (upsearch,downsearch)
#         regionString = '%s:%s:%s-%s' % (chrom,sense,start,stop)
#         print('plotting %s with window %s in region %s to outputfolder %s' % (geneName,plotName,regionString,outputFolder))
#         callGenePlot(dataFile,refID,plotName,annotFile,namesList,outputFolder,regionString,yScale)




#==========================================================================
#===================MAKE HYPER TABLE=======================================
#==========================================================================

# def hyperOccupancyEnrichment(dataFile,annotFile,namesList,tssGFFFile,transcribedGenesList,enrichedFolder,outputFolder,window=100000,macs=True):

#     dataDict = loadDataTable(dataFile)
        
#     startDict = makeStartDict(annotFile)
#     #making sure output folder exists
#     formatFolder(outputFolder,True)

#     #making a transcribed tss collection
#     tssGFF = parseTable(tssGFFFile,'\t')
    
#     #transcribedGenes
#     print('Loading transcribed genes')
#     transcribedGenesTable = parseTable(transcribedGenesList,'\t')
    
#     transcribedGenes = [int(line[0])-1 for line in transcribedGenesTable]

#     tssLoci = []
#     for i in transcribedGenes:
#         line = tssGFF[i]
#         tssLocus = Locus(line[0],line[3],line[4],line[6],line[1])
#         tssLoci.append(tssLocus)
#     print('making transcribed TSS collection')
#     tssCollection = LocusCollection(tssLoci,500)

#     for name in namesList:
#         print('Finding hyperoccupancy for %s' %(name))
#         hyperTable = []

#         if name.count('MM1S') == 1:
#             #mycTransLocusSense = Locus('chr14',105130609,105404388,'+','NM_002467')
#             #tssCollection.append(mycTransLocusSense)
#             mycTransLocusAnti = Locus('chr14',105404388-5000,105404388+5000,'-','NM_002467')
#             mycTransLocusSense = Locus('chr14',105130609-5000,105130609+5000,'+','NM_002467')

#             tssCollection.append(mycTransLocusAnti)
#             tssCollection.append(mycTransLocusSense)

    
#         peakTable = parseTable(enrichedFolder+dataDict[name]['enrichedMacs'],'\t')
#         ticker = 0
#         for peakLine in peakTable:
#             if ticker%1000 == 0:
#                 print(ticker)
#             ticker+=1
#             peakLocus = Locus(peakLine[0],peakLine[1],peakLine[2],'.',peakLine[3])
#             #check tss status
#             signal = float(peakLine[4])
#             proxGenes = []
#             overlappingTSSLoci = tssCollection.getOverlap(peakLocus,'both')
#             if len(overlappingTSSLoci) > 0:
#                 #this is a tss peak
#                 tss_peak = 1
#                 proxGenes+= [locus.ID() for locus in overlappingTSSLoci]
#             else:
#                 #this is an enhancer or outside peak
#                 tss_peak = 0

#                 peakCenter = (int(peakLine[1]) + int(peakLine[2]))/2
#                 searchLocus = Locus(peakLine[0],peakCenter-window,peakCenter+window,'.',peakLine[3])
#                 overlappingProxLoci = tssCollection.getOverlap(searchLocus,'both')
#                 proxGenes += [locus.ID() for locus in overlappingProxLoci]
#                 if peakLine[3] == 'MACS_peak_10552':
#                     print(overlappingProxLoci)
#                     print([locus.start() for locus in overlappingProxLoci])
#                     print([locus.ID() for locus in overlappingProxLoci])
#             if len(proxGenes) > 0:
#                 proxGenes = uniquify(proxGenes)

#                 proxString = string.join(proxGenes,',')
#                 proxNames = uniquify([startDict[geneID]['name'] for geneID in proxGenes])
#                 proxNamesString = string.join(proxNames,',')
#                 newLine = [peakLine[0],peakLine[1],peakLine[2],peakLine[3],signal,tss_peak,proxString,proxNamesString]
#                 hyperTable.append(newLine)

#         #now sort the hyper table
#         sortedHyperTable = [['CHROM','START','STOP','NAME','SIGNAL','TSS_PEAK','PROXIMAL_GENE_IDS','PROXIMAL_GENE_NAMES']]

#         peakOrder = order([float(line[4]) for line in hyperTable],decreasing=True)

#         for i in peakOrder:
#             sortedHyperTable.append(hyperTable[i])

#         #now do the gene assignment way

#         geneDict = {'totalSignal':defaultdict(float),'proximalEvents':defaultdict(int)}
#         for line in sortedHyperTable[1:]:
#             proximalGenes = line[6].split(',')
#             if len(proximalGenes) == 0:
#                 continue
#             else:
#                 signal = float(line[4])

#                 for geneID in proximalGenes:
#                     geneDict['proximalEvents'][geneID]+=1
#                     geneDict['totalSignal'][geneID]+=signal

#         geneCentricTable = [] 
#         geneList = geneDict['totalSignal'].keys()
#         print(len(geneDict['totalSignal'].keys()))
#         print(len(geneDict['proximalEvents'].keys()))
#         peakOrder = order([geneDict['totalSignal'][geneID] for geneID in geneList],decreasing=True)
        
#         for i in peakOrder:
#             geneID = geneList[i]
#             geneName = startDict[geneID]['name']
#             geneCentricTable.append([geneID,geneName,geneDict['proximalEvents'][geneID],geneDict['totalSignal'][geneID]])
        
#         collapsedGeneCentricTable = [['GENE_ID','GENE_NAME','PROXIMAL_EVENTS','SIGNAL']] 
        
#         usedNames = []
#         for line in geneCentricTable:
#             geneName = line[1]
#             if usedNames.count(geneName) == 0:
#                 collapsedGeneCentricTable.append(line)
#                 usedNames.append(geneName)
        
#         unParseTable(sortedHyperTable,'%s%s_hyperPeaks.txt' % (outputFolder,name),'\t')
#         unParseTable(collapsedGeneCentricTable,'%s%s_hyperGenes.txt' % (outputFolder,name),'\t')


#==========================================================================
#===================MAKE HYPER DENISTY=====================================
#==========================================================================
                
# def hyperOccupancyDensity(annotFile,hyperPeakFile,namesList,gffName,mappedFolder,outFolder,referenceName=''):

#     '''
#     finds the gene assignment of peaks
#     then maps to genes for each datasets in namesList
#     '''

#     startDict = makeStartDict(annotFile)
#     assignDict = {'promoter':defaultdict(list),'enhancer':defaultdict(list)}

#     hyperPeaks = parseTable(hyperPeakFile,'\t')
#     for line in hyperPeaks[1:]:
#         peakID = int(line[3].split('_')[-1])-1
#         geneList = line[6].split(',')
#         if int(line[5]) == 1:
#             assignDict['promoter'][peakID] += geneList
#         else:
#             assignDict['enhancer'][peakID] += geneList

#     if len(referenceName) == 0:
#         referenceName = namesList[0]

#     for name in namesList:
#         print('working on %s' % (name))
#         newTable = [['REFSEQ_ID','NAME','PROMOTER_REGIONS','PROMOTER_SIGNAL','DISTAL_REGIONS','DISTAL_SIGNAL']]
#         mappedGFF = parseTable('%s%s/%s_%s.gff' % (mappedFolder,gffName,gffName,name),'\t')
#         geneDict = {'promoter':defaultdict(list),'enhancer':defaultdict(list)}
#         for line in mappedGFF[1:]:
#             peakID = int(line[0].split('_')[-1])
#             regionSize = int(line[1].split(':')[-1].split('-')[1])-int(line[1].split(':')[-1].split('-')[0])
#             if assignDict['promoter'].has_key(peakID):
#                 for geneID in assignDict['promoter'][peakID]:
#                     geneDict['promoter'][geneID].append(float(line[2])*regionSize)
#             elif assignDict['enhancer'].has_key(peakID):
#                 for geneID in assignDict['enhancer'][peakID]:
#                     geneDict['enhancer'][geneID].append(float(line[2])*regionSize)
#             else:
#                 continue

        
#         if name == referenceName:
#             keeperIDs = []
#             allGenes = uniquify(geneDict['promoter'].keys()+geneDict['enhancer'].keys())
#             peakOrder = order([(sum(geneDict['enhancer'][x])+ sum(geneDict['promoter'][x])) for x in allGenes],decreasing=True)

#             peakOrderedGenes = [allGenes[x] for x in peakOrder]
            
#             usedNames =[]
#             for geneID in peakOrderedGenes:
#                 geneName = startDict[geneID]['name']
#                 if usedNames.count(geneName) == 1:
#                     continue
#                 usedNames.append(geneName)
#                 keeperIDs.append(geneID)

        
#         #now write the table!
#         for geneID in keeperIDs:
#             geneName = startDict[geneID]['name']
#             promoterRegions = len(geneDict['promoter'][geneID])
#             promoterSignal = sum(geneDict['promoter'][geneID])
#             enhancerRegions = len(geneDict['enhancer'][geneID])
#             enhancerSignal = sum(geneDict['enhancer'][geneID])


#             newLine = [geneID,geneName,promoterRegions,promoterSignal,enhancerRegions,enhancerSignal]
#             newTable.append(newLine)
#         print('writing table for %s' % (name))
#         unParseTable(newTable,'%s%s_%s_REF_geneMapped.txt' % (outFolder,name,referenceName),'\t')


# #making all of the ylfs

# #==========================================================================
# #===================CONVERTING SAMS TO YLFS================================
# #==========================================================================


# def makeYLFs(dataList=[],overwrite = False):
#     '''
#     makes sams for the dataset names specified. if blank, calls everything
#     '''

#     dataDict = loadDataTable()

#     if len(dataList) == 0:
#         dataList = dataDict.keys()
    
#     for name in dataList:
        
#         try:
#             sam = open(dataDict[name]['sam'],'r')
#         except IOError:
#             print('WARNING: No .sam file for %s' % (name))
#             continue
#         cmd = ' python /nfs/young_ata/CYL_code/samToYLF.py %s %s BWT1,BWT2' % (dataDict[name]['sam'],dataDict[name]['ylf'])
#         if overwrite:
#             print('making a .ylf for %s' % (name))
#             os.system(cmd)
#         else:
#             try:
#                 ylf = open(dataDict[name]['ylf'],'r')
#             except IOError:
#                 print('Making a .ylf for %s' % (name))
#                 os.system(cmd)



# #==========================================================================
# #===================CALLING THE ERROR MODEL================================
# #==========================================================================

        
# def callErrorModel(errorFolder,dataList = [],overwrite = False):

#     '''
#     for each dataset, calls the error model
#     '''

#     dataDict = loadDataTable()

#     if len(dataList) == 0:
#         dataList = dataDict.keys()

#     paramTemplate = [
#         ['TARGET READS FILE',''],
#         ['BACKGROUND READS FILE',''],
#         ['EXPERIMENT NAME', ''],
#         ['GENOME BUILD^(SPECIES AMD GENOME BUILD)', 'HG18'],
#         ['READ CATEGORIE(S) USED^(SEPARATE WITH COMMA)', 'BWT1,BWT2'],
#         ['MAXIMUM READ REPEATS', '2'],
#         ['GENOMIC BIN WIDTH^(BP)', '25'],
#         ['READ EXTENSION MODEL^(1=0to+200BP,2=-400BPto+600BP)', '1'],
#         ['P-VALUE THRESHOLD(S)^(SEPARATE WITH COMMA)', '1E-7,1E-8,1E-9'],
#         ['FRACTION OF GENOME AVAILABLE^(.5=50%)', '0.50'],
#         ['BIN TO REGION COMPRESS DISTANCE^(BP)', '200'],
#         ['COMPARISON FLANK DISTANCE^(BP)', '100'],
#         ['MINIMUM BIN ENRICHMENT OVER BACKGROUND^(NORMALIZED FOLD ENRICHMENT)', '2'],
#         ['MINIMUM REGION ENRICHMENT OVER BACKGROUND^(NORMALIZED FOLD ENRICHMENT)','5'],
#         ['MINIMUM REGION PEAK HEIGHT', '2'],
#         ['MINIMUM REGION LENGTH', '2'],
#         ['MODEL FOR CALLING GENES^(1=START SITES,2=FULL GENE)', '2'],
#         ['DISTANCE TO CALL GENES^(BP)', '2000'],
#         ['GENE LIST FILE(S)^(FILE LOCATION)','/nfs/young_ata/gene_tables/HG18_REFSEQ'],
#         ['WRITE BROWSER TRACK^(0=NO,1=YES,2=WRITE BACKGROUND TRACK ALSO)', '2'],
#         ['BROWSER TRACK FLOOR^(COUNTS)', '2'],
#         ['NORMALIZE TO READS PER MILLION^(1=YES,0=NO)', '0'],
#         ['DENSITY PLOT REGION(S)^(NAME,CHR,START,END)', '0'],
#         ['METAGENE PARAMETERS^(UPSTREAM_BP,DOWNSTREAM_BP,GENOMIC_BINS_PER_CLUSTERGRAM_BIN,NUMBER_OF_INGENE_BINS)','2000,2500,2,50'],
#         ['WRITE BIN COUNT HISTOGRAM^(1=YES,0=NO)', '1'],
#         ['WRITE DENSITY FILE^(1=YES,0=NO)', '0'],
#         ['WRITE GENOME SEQUENCE AT ENRICHMENT PEAKS^(+/- BP WINDOW)', '0'],
#         ['OUTPUT FOLDER LOCATION^(0=CURRENT DIRECTORY)', '0'],
#         ['END']
#         ]
    
#     for name in dataList:
        
#         #don't run error model on background datasets
#         if dataDict[name]['background'] == 'NONE':
#             continue
        
#         paramFile = list(paramTemplate)
#         paramFile[0][1] = dataDict[name]['ylf']
#         backgroundName =  dataDict[name]['background']
#         paramFile[1][1] = dataDict[backgroundName]['ylf']
#         paramFile[2][1] = name
        
#         #if overwrite is false, check for an enriched region file and skip if it exists
#         if overwrite == False:

#             if dataDict[name]['enriched'] != 'NONE':
#                 continue

#         #write the parameter file
#         paramFileName = name + '_ChIPseq_totalcounts_enrichedregions.txt'
#         print(paramFileName)
#         unParseTable(paramFile,errorFolder+paramFileName,'\t')

#         #now call the error model for this file
        
#         print('changing directory to %s' % (errorFolder))
#         os.chdir(errorFolder)
#         cmd = "bsub -o /dev/null -R 'rusage[mem=2200]' 'python /nfs/young_ata/GMF_code/Solexa/ENRICHED_REGIONS_CODE_NEW/SOLEXA_ENRICHED_REGIONS.py %s'" % (paramFileName)
#         print(cmd)
#         os.system(cmd)


# #==========================================================================
# #===================FORMATTING ERROR MODEL=================================
# #==========================================================================

# def formatErrorModelOutput(dataTableFile,errorOutput,wiggleFolder,enrichedRegionsFolder):

#     '''
#     takes a folder of error model output and moves the enriched regions and wiggles to the right folder
#     updates the dataTable with the names of the enrichedRegions
#     '''

#     dataTable = parseTable(dataTableFile,'\t')
        
#     nameList = [line[3] for line in dataTable[1:]]

#     if errorOutput[-1] != '/':
#         errorOutput+='/'

#     if wiggleFolder[-1] != '/':
#         wiggleFolder+='/'

#     if enrichedRegionsFolder[-1] != '/':
#         enrichedRegionsFolder+='/'

#     outputFolders = os.listdir(errorOutput)

#     #filter for folders

#     outputFolders = filter(lambda x: x.split('.')[-1] != 'txt',outputFolders)

#     #get rid of any hidden files

#     outputFolders = filter(lambda x: x[0] != '.',outputFolders)

#     for output in outputFolders:

#         #get the name of the dataset
#         name = string.join(output.split('_')[:-2],'_')

#         #first move the wiggle
#         cmd = 'mv %s%s/*.WIG.gz %s &' % (errorOutput,output,wiggleFolder)
#         print(cmd)
#         os.system(cmd)
        
#         #next move the enriched region
#         cmd = 'mv %s%s/1e-09/ENRICHED_REGIONS_%s %s &' % (errorOutput,output,output,enrichedRegionsFolder)
#         print(cmd)
#         os.system(cmd)

#         #next update the dataTable
#         i = nameList.index(name)
#         dataTable[i+1][5] = 'ENRICHED_REGIONS_'+output
#         print(dataTable[i+1])

#     unParseTable(dataTable,dataTableFile,'\t')



# #==========================================================================
# #===================MAKING GFFS OF TSS REGIONS=============================
# #==========================================================================
        
# def makeTSSGFFs(dataTable,annotFile,upstream,downstream,gffFolder,species='HG18'):

#     '''
#     we're making several kinds of gffs for this study
#     all regions will be 10kb in size,
#     all binding sites will be within 2kb of one another or reference points
#     '''
#     dataDict = loadDataTable(dataTable)

#     if gffFolder[-1] != '/':
#         gffFolder+='/'



#     #TSS
#     #tss gff for all genes
#     tssGFF = makeTSSGFF(annotFile,upstream,downstream)
    
#     unParseTable(tssGFF,gffFolder + '%s_TSS_ALL_-%s_+%s.gff' % (species,upstream,downstream),'\t')

# #==========================================================================
# #===================MAKING GFFS OF GENE BODY REGIONS=======================
# #==========================================================================
        
# def makeBodyGFFs(annotFile,gffFolder):

#     '''
#     we're making several kinds of gffs for this study
#     all regions will be 10kb in size,
#     all binding sites will be within 2kb of one another or reference points
#     '''
#     dataDict = loadDataTable()

#     if gffFolder[-1] != '/':
#         gffFolder+='/'

#     startDict = makeStartDict(annotFile)

#     #TSS
#     #tss gff for all genes

#     bodyGFF = []

#     for gene in startDict.keys():

#         chrom = startDict[gene]['chr']
#         sense = startDict[gene]['sense']
#         start = startDict[gene]['start'][0]
#         end = startDict[gene]['end'][0]

#         if sense == '+':

#             bodyGFF.append([chrom,gene,'',start+300,end+3000,'',sense,'',gene])

#         else:

#             bodyGFF.append([chrom,gene,'',end-3000,start-300,'',sense,'',gene])
    
#     unParseTable(bodyGFF,gffFolder + 'HG18_BODY_ALL_+300_+3000.gff','\t')



# #==========================================================================
# #===================MAKING GFFS OF ENHANCER REGIONS========================
# #==========================================================================


# def makeEnhancerGFFs(dataTable,gffName,nameList,annotFile,gffFolder,enrichedFolder,macs=True):
#     '''
#     find all possible enhancers.
#     enhancers defined as h3k27ac binding sites +/-5kb outside of promoters
#     we define center of enhancer as the center of the bound region
#     '''

#     dataDict = loadDataTable(dataTable)


#     if enrichedFolder[-1] != '/':
#         enrichedFolder+='/'

#     if gffFolder[-1] != '/':
#         gffFolder+='/'

#     #nameList = ['H128_H3K27AC','H2171_H3K27AC','MM1S_H3K27AC_DMSO','MM1S_H3K27AC_JQ1','U87_H3K27AC','P493-6_T0_H3K27AC','P493-6_T1_H3K27AC','P493-6_T24_H3K27AC']

#     #first make the tss collection
#     tssGFF = makeTSSGFF(annotFile,5000,5000)
#     tssCollection = gffToLocusCollection(tssGFF)

#     #make a blank collection to load enhancers into
#     enhancerCollection = LocusCollection([],500)

#     #don't allow overlapping enhancers
#     for name in nameList:


#         print('finding enhancers in %s' % (name))
#         if macs:
#             boundCollection = importBoundRegion(enrichedFolder + dataDict[name]['enrichedMacs'],name)
#         else:
#             boundCollection = importBoundRegion(enrichedFolder + dataDict[name]['enriched'],name)
#         for locus in boundCollection.getLoci():
#             #make sure an overlapping enhancer doesn't already exist
#             if len(tssCollection.getOverlap(locus,'both')) == 0 and len(enhancerCollection.getOverlap(locus,'both')) == 0:
#                 center = (locus.start()+locus.end())/2
#                 gffLocus = Locus(locus.chr(),center-5000,center+5000,'.',locus.ID())
#                 enhancerCollection.append(gffLocus)


#     enhancerGFF = locusCollectionToGFF(enhancerCollection)
#     print('Found %s enhancers in %s' % (len(enhancerGFF),gffName))
#     unParseTable(enhancerGFF,gffFolder+'HG18_ENHANCERS_%s_-5000_+5000.gff' % (gffName),'\t')


# #==========================================================================
# #===================MAKE ELEMENTS GFF======================================
# #==========================================================================

# def makeElementsGFF(tssGFFFile,transcribedListFile,enhancerGFFFile,enhancerListFile,rnaTableFile,gffFolder,elementsName):


#     '''
#     makes a gff of elements active promoters, silent promoters, active enhancers, rRNA genes, tRNA genes
#     '''


#     elementGFF = []

#     #start with the promoters

#     transcribedList= parseTable(transcribedListFile,'\t')

#     #all lists are 1 indexed
#     transcribedList = [int(line[0])-1 for line in transcribedList]

#     tssGFF = parseTable(tssGFFFile,'\t')

#     for i in range(len(tssGFF)):

#         newLine = list(tssGFF[i])

#         if transcribedList.count(i) == 1:
#             newLine[1] = 'ACTIVE'
#             transcribedList.remove(i)
#         else:
#             newLine[1] = 'SILENT'
        
#         newLine[3] = int(newLine[3]) +4000
#         newLine[4] = int(newLine[4]) -4000

#         elementGFF.append(newLine)

#     #now lets get enhancers

#     enhancerList = parseTable(enhancerListFile,'\t')
    
    
#     enhancerList = [int(line[0])-1 for line in enhancerList]


#     enhancerGFF = parseTable(enhancerGFFFile,'\t')

#     for i in enhancerList:

#         newLine = list(enhancerGFF[i])
#         newLine[1] = 'ENHANCER'
#         newLine[3] = int(newLine[3]) +4000
#         newLine[4] = int(newLine[4]) -4000

#         elementGFF.append(newLine)


#     #now lets do RNA genes

#     rnaTable = parseTable(rnaTableFile,'\t')

#     for line in rnaTable[1:]:
        
#         #check if it's a tRNA
#         if line[6] == 'Eddy-tRNAscanSE' and float(line[8]) >20:
#             #now we know we have a potential tRNA
#             #check the score to make sure it's legit
#             #use the 20 cutoff

#             if line[5] == '+':

#                 newLine = [line[0],'tRNA','',int(line[1]) - 1000,int(line[1]) + 1000,float(line[8]),'+','',line[3]]
#             else:
#                 newLine = [line[0],'tRNA','',int(line[2]) - 1000,int(line[2]) + 1000,float(line[8]),'+','',line[3]]                    
#             elementGFF.append(newLine)
#         #now lets get rRNA genes
#         if line[6] == 'Eddy-BLAST-otherrnalib' and line[7] == 'rRNA':
#             if line[5] == '+':

#                 newLine = [line[0],'rRNA','',int(line[1]) - 1000,int(line[1]) + 1000,float(line[8]),'+','',line[3]]
#             else:
#                 newLine = [line[0],'rRNA','',int(line[2]) - 1000,int(line[2]) + 1000,float(line[8]),'+','',line[3]]
#             elementGFF.append(newLine)
#     unParseTable(elementGFF,gffFolder+elementsName,'\t')

# #==========================================================================
# #===================MAKING GFF FROM MACS PEAKS=============================
# #==========================================================================


# def makePeakGFF(namesList,macsFolder,window,gffFolder):

#     '''
#     takes a macs summit file and makes a gff of all summits +/- window and spits it out in the gff folder
#     '''

#     dataDict = loadDataTable()
#     if macsFolder[-1] != '/':
#         macsFolder+='/'

#     if gffFolder[-1] !='/':
#         gffFolder+='/'

#     for name in namesList:

#         summitBed = parseTable(macsFolder+name+'/' +name+'_summits.bed','\t')

#         summitGFF = []

#         for line in summitBed:

#             chrom = line[0]
#             start = int(line[1]) - window
#             end = int(line[1]) + window
#             name = line[3]
#             score = line[4]

#             newLine = [chrom,name,'',start,end,score,'.','',name]

#             summitGFF.append(newLine)

#         unParseTable(summitGFF,'%sHG18_%s_summits_-%s_+%s.gff' % (gffFolder,name,window,window),'\t')



# #==========================================================================
# #===================RANKING E-BOXESFROM MACS PEAKS=========================
# #==========================================================================


# def rankEboxes(name,genomeDirectory,window,macsFolder,motifFolder,asDict=False):

#     '''
#     takes a summit bed and ranks eboxes by height from a sequence in a window around the bed

#     '''
    
#     #hardcoded ebox motif

#     motifRegex = 'CA..TG'
#     if macsFolder[-1] != '/':
#         macsFolder+='/'

#     if motifFolder[-1] != '/':
#         motifFolder+='/'

#     dataDict = loadDataTable()
#     window = int(window)

#     summitBed = parseTable(macsFolder+name+'/' +name+'_summits.bed','\t')
#     eboxDict = defaultdict(list)
#     ticker= 0
#     for line in summitBed:
#         if ticker % 1000 == 0:
#             print(ticker)
#         ticker+=1

#         chrom = line[0]
#         peakName = line[3]
#         sense = '.'


#         start = int(line[1])-window
#         end = int(line[1])+window
#         height = float(line[4])

#         sequenceLine = fetchSeq(genomeDirectory,chrom,start,end,True)
        
#         motifVector = []
#         matches = re.finditer(motifRegex,string.upper(sequenceLine))
#         if matches:
#             for match in matches:
#                 motifVector.append(match.group())
        
#         #count only 1 of each motif type per line
#         motifVector = uniquify(motifVector)
#         for motif in motifVector:

#             eboxDict[motif].append(height)


#     eboxTable =[]
#     eboxTableOrdered =[['EBOX','OCCURENCES','AVG_HEIGHT']]
#     for ebox in eboxDict.keys():
#         newLine = [ebox,len(eboxDict[ebox]),mean(eboxDict[ebox])]
#         eboxTable.append(newLine)


#     occurenceOrder = order([line[2] for line in eboxTable],decreasing=True)
    
#     for x in occurenceOrder:
#         eboxTableOrdered.append(eboxTable[x])
#     print(eboxTableOrdered)
#     unParseTable(eboxTableOrdered,motifFolder+name+'_eboxes_unique.txt','\t')


# #==========================================================================
# #===================MAKE A FASTA FROM A PEAK FILE==========================
# #==========================================================================
    

        
# def makePeakFastas(name,genomeDirectory,window,tssGFFFile,enhancerGFFFile,macsFolder,outFolder,top=''):

#     '''
#     makes a ranked fasta based on peak height. best = True takes only the best peak from each region
#     '''

#     if macsFolder[-1] != '/':
#         macsFolder+='/'

#     if outFolder[-1] != '/':
#         outFolder+='/'

#     dataDict = loadDataTable()
#     window = int(window)

        
#     #load in the tss GFF and enhancerGFF
#     tssGFF = parseTable(tssGFFFile,'\t')
#     tssCollection = gffToLocusCollection(tssGFF)

#     enhancerGFF = parseTable(enhancerGFFFile,'\t')
#     enhancerCollection = gffToLocusCollection(enhancerGFF)

#     peakTable = parseTable(macsFolder+name+'/' +name+'_summits.bed','\t')

#     tssFastaList = []
#     enhancerFastaList = []
#     ticker =0
#     outsideBoth = 0

#     for line in peakTable:
#         if ticker % 1000 == 0:
#             print(ticker)
#         ticker+=1

#         chrom = line[0]
#         peakName = line[3]
#         sense = '.'


#         start = int(line[1])-window
#         end = int(line[1])+window
#         height = float(line[4])

#         peakLocus = Locus(chrom,start,end,sense,peakName)
#         #check to see if the peak is near a promoter or enhancer
#         if len(tssCollection.getOverlap(peakLocus,'both')) > 0 and len(enhancerCollection.getOverlap(peakLocus,'both')) == 0:
        
#             fastaHeader = '>%s:%s:%s-%s|%s|%s' % (chrom,sense,start,end,peakName,height)
#             fastaLine = fetchSeq(genomeDirectory,chrom,start,end,True)
#             tssFastaList.append([fastaHeader,fastaLine])
#         if len(tssCollection.getOverlap(peakLocus,'both')) == 0 and len(enhancerCollection.getOverlap(peakLocus,'both')) >0:

        
#             fastaHeader = '>%s:%s:%s-%s|%s|%s' % (chrom,sense,start,end,peakName,height)
#             fastaLine = fetchSeq(genomeDirectory,chrom,start,end,True)
#             enhancerFastaList.append([fastaHeader,fastaLine])

#         if len(tssCollection.getOverlap(peakLocus,'both')) == 0 and len(enhancerCollection.getOverlap(peakLocus,'both')) == 0:
#             outsideBoth+=1


#     #this returns a fastaList
#     print('found %s tss bound regions' % len(tssFastaList))
#     print('found %s enhancer bound regions' % len(enhancerFastaList))
#     print('found %s bound regions outside both' % outsideBoth)
    
    
#     #write the tss fasta
#     peakOrder = order([float(fasta[0].split('|')[-1]) for fasta in tssFastaList],decreasing=True)

#     outFasta = open('%s%s_TSS_-%s_+%s_top%s.fasta' % (outFolder,name,window,window,top),'w')
#     if not top:
#         top = len(tssFastaList)
#     else:
#         top = int(top)

#     for i in range(top):

#         peakIndex = peakOrder[i]

#         outFasta.write(tssFastaList[peakIndex][0]+'\n')
#         outFasta.write(tssFastaList[peakIndex][1]+'\n')


#     outFasta.close()

#     #write the enhancer fasta

#     if not top:
#         top = len(enhancerFastaList)
#     else:
#         top = int(top)

#     peakOrder = order([float(fasta[0].split('|')[-1]) for fasta in enhancerFastaList],decreasing=True)

#     outFasta = open('%s%s_ENHANCER_-%s_+%s_top%s.fasta' % (outFolder,name,window,window,top),'w')

#     for i in range(top):

#         peakIndex = peakOrder[i]

#         outFasta.write(enhancerFastaList[peakIndex][0]+'\n')
#         outFasta.write(enhancerFastaList[peakIndex][1]+'\n')


#     outFasta.close()            

# #==========================================================================
# #===================RANKING E-BOXES FROM FASTA=============================
# #==========================================================================

# def rankEboxFasta(fastaFile,eboxTableFile,motifFolder):

#     '''
#     finds all of the eboxes in a fasta and writes a table of how often they occur and the avg. score
#     for that particular sequence
#     '''

#     motifRegex = 'CA..TG'
#     if motifFolder[-1] != '/':
#         motifFolder+='/'

#     fastaName = fastaFile.split('/')[-1].split('.')[0]

#     fasta = parseTable(fastaFile,'\t')


#     eboxDict = defaultdict(float)

#     eboxTable = parseTable(eboxTableFile,'\t')


#     for line in eboxTable[1:]:
#         eboxDict[line[0]] = float(line[2])
#     print(eboxDict)

#     occurenceDict = defaultdict(int)

#     for line in fasta:
#         line = line[0]
#         if line[0] == '>':
#             continue


#         matches = re.finditer(motifRegex,string.upper(line))
#         if matches:
#             for match in matches:
#                 occurenceDict[match.group()]+=1


#     strengthOrder = order([line[2] for line in eboxTable[1:]],decreasing=True)
#     print(strengthOrder)
#     eboxList = [eboxTable[x+1][0] for x in strengthOrder]

#     print(occurenceDict)
#     occurenceTable = [['EBOX','OCCURENCES']]

#     for ebox in eboxList:

#         occurenceTable.append([ebox,occurenceDict[ebox]])


#     unParseTable(occurenceTable,motifFolder+fastaName+'_eboxes.txt','\t')


    





                                                                                                            
# #==========================================================================
# #===================CALL MEME==============================================
# #==========================================================================

# def callMEME(fastaFolder,overwrite = False):

#     '''
#     calls MEME on a folder full of fastas. has option to overwrite
#     '''
#     if fastaFolder[-1]!='/':
#         fastaFolder+='/'
#     fastaFolderList = os.listdir(fastaFolder)

#     fastaFileList = filter(lambda x: x.split('.')[-1] == 'fasta',fastaFolderList)

#     for fastaFile in fastaFileList:

#         fastaName = fastaFile.split('.')[0]
#         #the don't overwrite condition
#         if overwrite == False and fastaFolderList.count(fastaName) == 1:
            
#             print('yay')
#             continue
#         print('calling MEME on %s' % fastaName)
#         output = fastaFolder+fastaName+'/'

#         cmd = "bsub -R 'rusage[mem=2200]' 'meme -dna -evt 1 -mod zoops -nmotifs 10 -minw 4 -maxw 10 -revcomp -o %s -maxsize 1000000 %s'" % (output,fastaFolder+fastaFile)
#         print(cmd)
#         os.system(cmd)

    
    
# #==========================================================================
# #===================MAPPING BAMS TO GFFS===================================
# #==========================================================================

# def mapMycBams(cellTypeList,gffList,mappedFolder,dataFile = '',nBin = 200,overWrite =False,nameList = []):
    
#     '''
#     for each gff maps all of the data and writes to a specific folder named after the gff
#     '''
#     dataDict = loadDataTable(dataFile)
#     if mappedFolder[-1] != '/':
#         mappedFolder+='/'

#     for gffFile in gffList:
        
#         #check to make sure gff exists
#         try:
#             foo = open(gffFile,'r')
#         except IOError:
#             print('ERROR: GFF FILE %s DOES NOT EXIST' % (gffFile))

#         gffName = gffFile.split('/')[-1].split('.')[0]
        
#         #make this directory
#         try:
#             foo = os.listdir(mappedFolder+gffName)
#         except OSError:
#             os.system('mkdir %s%s' % (mappedFolder,gffName))

#         outdir = mappedFolder+gffName+'/'

#         if len(nameList) == 0:
#             nameList = dataDict.keys()
        

#         for name in nameList:
            
#             #filter based on celltype
#             cellName = name.split('_')[0]
#             if cellTypeList.count(cellName) != 1:
#                 continue
#             fullBamFile = dataDict[name]['bam']
#             outFile = outdir+gffName+'_'+name+'.gff'

#             if overWrite:
#                 cmd1 = "bsub -o /dev/null -R 'rusage[mem=2200]' 'python /nfs/young_ata/scripts/bamToGFF.py -u -d -f 1 -e 200 -r -m %s -b %s -i %s -o %s'" % (nBin,fullBamFile,gffFile,outFile)
#                 print(cmd1)
#                 os.system(cmd1)
#             else:
#                 try:
#                     Foo = open(outFile,'r')
#                     print('File %s Already Exists, not mapping' % (outFile))
#                 except IOError:

#                     cmd1 = "bsub -o /dev/null -R 'rusage[mem=2200]' 'python /nfs/young_ata/scripts/bamToGFF.py -u -d -f 1 -e 200 -r -m %s -b %s -i %s -o %s'" % (nBin,fullBamFile,gffFile,outFile)
#                     print(cmd1)
#                     os.system(cmd1)

# #==========================================================================
# #===================PULLING MOTIFS FROM FASTAS=============================
# #==========================================================================


# def rankedFastaMotifs(fastaFile,outFolder):

    
#     '''
#     takes a ranked fasta and spits out a bunch of files that are nifty
#     '''

#     if outFolder[-1] != '/':
#         outFolder+='/'
#     fastaName = fastaFile.split('/')[-1].split('.')[0]
#     try:
#         os.listdir('%s%s_rankedMotifs/' % (outFolder,fastaName))
#     except  OSError:
#         os.system('mkdir %s%s_rankedMotifs/' % (outFolder,fastaName))
#     fasta = parseTable(fastaFile,'\t')
#     #hard coded for E-boxes
#     motif = '[ACTG]CA[ACTG]{2}TG[ATCG]'
#     motifLen = 8
#     nBins = 10
    
#     binSize = len(fasta)/nBins

#     heightTable = [['BIN','HEIGHT','NSITES','A','C','G','T']]
    

#     for i in range(nBins):
#         compDict = defaultdict(int)
#         start = i*binSize
#         end = (i+1) * binSize
#         heightVector = []
#         motifVector = []
#         for j in range(start,end,2):
#             headerLine = fasta[j][0]
#             sequenceLine = fasta[j+1][0]
#             heightVector.append(float(headerLine.split('|')[-1]))
#             matches = re.finditer(motif,string.upper(sequenceLine))
#             if matches:
#                 for match in matches:
#                     motifVector.append(match.group())
#             lineComp = composition(sequenceLine,['A','T','C','G'])
#             for x in lineComp.keys():
#                 compDict[x] += lineComp[x]
#         totalSequence = sum([compDict[x] for x in compDict.keys()])
        
#         compLine = [float(compDict[x])/totalSequence for x in ['A','C','G','T']]

#         heightTable.append([i+1,mean(heightVector),len(motifVector)]+compLine)
#         motifMatrix = [['A','C','G','T']]
#         for position in range(motifLen):
#             positionVector = [x[position] for x in motifVector]
#             positionCounts = [positionVector.count('A'),positionVector.count('C'),positionVector.count('G'),positionVector.count('T')]
#             positionCounts = [round(100*float(x)/len(positionVector),0) for x in positionCounts]
#             motifMatrix.append(positionCounts)
#         unParseTable(motifMatrix,'%s%s_rankedMotifs/motif_nBin_%s.txt' % (outFolder,fastaName,i),'\t')
#         motifFasta = [[x] for x in motifVector]
#         unParseTable(motifFasta,'%s%s_rankedMotifs/motif_nBin_%s.fasta' % (outFolder,fastaName,i),'\t')



#     unParseTable(heightTable,'%s%s_rankedMotifs/motif_heights.txt' % (outFolder,fastaName),'\t')

# #==========================================================================
# #===================GFF MOTIF TO BED ======================================
# #==========================================================================

# def gffMotifToBed(name,motif,genomeDirectory,gffFile,bedFolder):

#     '''
#     writes a bed track showing where motifs are within gff regions
#     '''

#     gffName = gffFile.split('/')[-1].split('.')[0]
#     trackLine = 'track name="%s" description="%s" visibility=3' % (name,gffName)

#     bed = open('%s%s_%s.bed' % (bedFolder,gffName,name),'w')
#     bed.write(trackLine+'\n')
#     gff = parseTable(gffFile,'\t')
#     ticker =0
#     for line in gff:
#         if ticker%1000 == 0:
#             print ticker
#         ticker+=1
#         chrom = line[0]
#         start = int(line[3])
#         end = int(line[4])
#         try:
#             sequence = fetchSeq(genomeDirectory,chrom,start,end,True)
#         except IOError:
#             continue
#         motifs = re.finditer(motif,string.upper(sequence))
#         for match in motifs:
    
#             newLine = string.join([chrom,str(start+match.start()),str(start+match.end())],'\t')
#             bed.write(newLine+'\n')


#     bed.close()

    
# #check for mapped bams
# #==========================================================================
# #===================CHECK MAPPED BAMS======================================
# #==========================================================================

# def checkMappedBams(gffList,mappedFolder):

#     '''
#     goes through the mapped folder to make sure each bam was mapped to each gff correctly
#     prints an error message if things messed up
#     '''
#     dataDict = loadDataTable()
#     if mappedFolder[-1] != '/':
#         mappedFolder+='/'

#     for gffFile in gffList:

#         gffName = gffFile.split('/')[-1].split('.')[0]
        
#         outdir = mappedFolder+gffName+'/'

#         nameList = dataDict.keys()
        
#         for name in nameList:

#             outFile = outdir+gffName+'_'+name+'.gff'
#             try:
#                 foo = open(outFile,'\r')
#             except IOError:
#                 print('NO MAPPED BAM FOR %s FOUND FOR DATASET %s' % (gffName,name))



        

# #==========================================================================
# #===================FINDING TARGET GENES===================================
# #==========================================================================



   
# def mergeTargetLists(startDict,namesList,enrichedFolder):
#     dataDict = loadDataTable()
#     targetGeneList = []


#     for name in namesList:

#         enrichedFile = enrichedFolder + dataDict[name]['enriched']
#         targetGeneList +=  targetGenes(startDict,enrichedFile,-1000,1000)

#     targetGeneList = uniquify(targetGeneList)
#     return targetGeneList


# #==========================================================================
# #===================MAP ENRICHED REGIONS TO GFF============================
# #==========================================================================

# def mapEnrichedToGFF(dataFile,setName,gffList,cellTypeList,enrichedFolder,mappedFolder,macs=True):

#     '''
#     maps enriched regions from a set of cell types to a set of gffs
#     tries to make a new folder for each gff
#     '''

#     dataDict = loadDataTable(dataFile)
#     if enrichedFolder[-1] != '/':
#         enrichedFolder+='/'

#     if mappedFolder[-1] != '/':
#         mappedFolder+='/'

#     for gffFile in gffList:

#         gffName = gffFile.split('/')[-1].split('.')[0]
#         print('making enriched regions to %s' % (gffName))
#         #make the gff into a collection


        
#         outdir = mappedFolder+gffName+'/'

#         try:
#             foo = os.listdir(mappedFolder+gffName)
#         except OSError:
#             os.system('mkdir %s' % (outdir))
        
#         #first filter the name list
#         cellTypeNameList =[] 

#         for name in dataDict.keys():

#             #check to make sure in the right celltype
#             #also make sure to not process WCEs
#             if dataDict[name]['background'] == 'NONE':
#                 continue
#             cellName = name.split('_')[0]
#             if macs == True:
#                 if cellTypeList.count(cellName) == 1 and dataDict[name]['enrichedMacs'] != 'NONE':
#                     cellTypeNameList.append(name)

#             else:
#                 if cellTypeList.count(cellName) == 1 and dataDict[name]['enriched'] != 'NONE':
#                     cellTypeNameList.append(name)

#         cellTypeNameList.sort()

#         mappedGFF = [['GFF_LINE','ID'] + cellTypeNameList]
#         #now we go through the gff and fill in stuff
#         gffTable = parseTable(gffFile,'\t')

#         gffLoci = []
#         for line in gffTable:
#             gffLocus = Locus(line[0],line[3],line[4],line[6],line[8])
#             gffLine = gffLocus.__str__()
#             gffID = line[1]
            
#             gffLoci.append(gffLocus)
#             mappedGFF.append([gffLine,gffID])
            
#         for name in cellTypeNameList:
#             print('dataset %s' % (name))
#             if macs:
#                 enrichedCollection = importBoundRegion(enrichedFolder + dataDict[name]['enrichedMacs'],name)
#             else:
#                 enrichedCollection = importBoundRegion(enrichedFolder + dataDict[name]['enriched'],name)
#             for i in range(len(gffLoci)):
#                 if len(enrichedCollection.getOverlap(gffLoci[i],'both')) > 0:
#                     mappedGFF[i+1].append(1)
#                 else:
#                     mappedGFF[i+1].append(0)


#         unParseTable(mappedGFF,outdir+gffName+'_'+setName+'.txt','\t')


# #==========================================================================
# #===================FORMATTING NANOSTRING DATA=============================
# #==========================================================================

# def formatNanoString(nanoStringFolder,sampleNames,tssGFFFile,geneListFile,output):

#     '''
#     takes a folder of nanostring data numbered in our current schema
#     and formats it into a pickle and a table
#     '''

#     #first get the tss gff file loaded

#     tssGFF = parseTable(tssGFFFile,'\t')

#     geneListTable = parseTable(geneListFile,'\t')
    
#     geneList = [int(line[0]) -1 for line in geneListTable]

#     activeList = [tssGFF[i][1] for i in geneList]



#     nanoStringFileList = os.listdir(nanoStringFolder)

#     nanoDict = {'0hr':defaultdict(list),
#                 '1hr':defaultdict(list),
#                 '24hr':defaultdict(list),
#                 'notet':defaultdict(list),
#                 'H128':defaultdict(list),
#                 'H2171':defaultdict(list),
#                 }

#     sampleDict = {'01':'0hr',
#                   '02':'0hr',
#                   '03':'1hr',
#                   '04':'1hr',
#                   '05':'24hr',
#                   '06':'24hr',
#                   '07':'notet',
#                   '08':'notet',
#                   '09':'H2171',
#                   '10':'H2171',
#                   '11':'H128',
#                   '12':'H128',
#                   }

#     nameDict = defaultdict(str)
#     activityDict = defaultdict(str)
#     for nanoFile in nanoStringFileList:
#         if nanoFile.split('.')[-1] != 'RCC':
#             continue
#         sampleID = nanoFile.split('.')[0].split('_')[-1]
#         sampleName = sampleDict[sampleID]

#         dataTable = parseTable(nanoStringFolder+nanoFile,'\r')

#         for line in dataTable:

#             if line[0].count('NM') == 1:
#                 geneName = line[0].split(',')[1]
#                 refseqID = line[0].split(',')[2].split('.')[0]
#                 nameDict[refseqID] = geneName
#                 if activeList.count(refseqID) == 1:
#                     activityDict[refseqID] ='ACTIVE'
#                 else:
#                     activityDict[refseqID] = 'SILENT'
#                 geneCounts = int(line[0].split(',')[-1])
#                 nanoDict[sampleName][refseqID].append(geneCounts)


#     nanoTable = [['REFSEQ_ID','GENE','ACTIVITY'] +sampleNames]
#     nanoGeneList = nanoDict['0hr'].keys()


#     for refseqID in nanoGeneList:

#         newLine = [refseqID,nameDict[refseqID],activityDict[refseqID]]+[sum(nanoDict[x][refseqID][0:2]) for x in sampleNames]
#         nanoTable.append(newLine)



#     unParseTable(nanoTable,output,'\t') 

   




# #==========================================================================
# #===================PLOTTING INDIVIDUAL GENES=============================
# #==========================================================================

# def makeGeneListFile(gffBindingFile,setList,output):

#     '''

#     AND operators within lists, OR operators outside of lists
#     [[A,B],[C,D]] = (A AND B) OR (C AND D) for this row
#     [[A],[B],[C],[D]] = A OR B OR C OR D for this row

#     '''
#     geneListFile = []
#     boundGFFTable = parseTable(gffBindingFile,'\t')
#     header = boundGFFTable[0]
    
#     #convert the setList into column numbers
#     columnSet = []
#     for bindingSet in setList:
#         try:
#             columnSet.append([header.index(x) for x in bindingSet])
#         except ValueError:
#             print('ERROR: not all datasets in binding table')
#             exit()

#     for i in range(1,len(boundGFFTable),1):
        
#         line = boundGFFTable[i]
        
#         #print(columnSet)
#         #print(i)
#         #if any of these end up being true, the line gets added
#         for andColumns in columnSet:
            
#             bindingVector = [int(line[x]) for x in andColumns]
#             #print(bindingVector)
#             if bindingVector.count(1) == len(bindingVector):
#                 geneListFile.append(i)
#                 break
#     print(len(geneListFile))
#     unParseTable(geneListFile,output,'')

#==========================================================================
#===================MAKE UCSC TRACKHUB FILES===============================
#==========================================================================

def makeTrackHub(analysis_name,project_folder,chrom_sizes, dataFileList=[], wiggle_dir='',web_dir='/storage/cylin/web/Lin_Lab_Track_Hubs/',hub_name='',hub_short_lab='',hub_long_lab='',EMAIL='',fileType='bigWig',col='0,0,0',scaled=False):
    #dataFileList will take several data tables and use them to create a track hub. This will not include background samples, because they do not have wiggle files
    #analysis_name is the name of your project
    #wiggle_dir is the path to where the wiggle files for this analysis live; will default to '/storage/cylin/grail/projects/analysis_name/wiggles' if left blank
    #chrom_sizes is a required file for wigToBigWig. You can make these by running fetchChromSizes on your assembly. These are also available for download from UCSC, too.
    #web_dir is the path to the outward facing directory where your data will be accessible to the genome browser
    #hub_name will default to the analysis name if left empty
    #hub_short_lab is a short label description of the hub and defaults to the analysis name if left blank
    #hub_long_lab is a long label description of the hub and will also default to the analysis name if left blank
    #EMAIL is your e-mail address
    #fileType is set to bigWig by default for this script, but trackhub supports several other file types. This can be updated in the future.
    #col is the color of the tracks; default black. This can be edited on the UCSC genome browser as well


    folderName = '{}{}'.format(
    web_dir,
    analysis_name
    )

    formatFolder(folderName,create=True)

    group = web_dir.split('/')[2]
    track_folder = web_dir.split('/web')[1]

    urlBase = 'http://taco-wiki.grid.bcm.edu/'
    if wiggle_dir == '':
        wiggle_dir = project_folder+'wiggles/'

    allData = []
    bgNames = []

    #creates a large list from all data files selected to create track hub; removes background samples
    for file in dataFileList:
        dataTable = utils.parseTable(file,'\t')
        for line in dataTable[1:]:
            allData.append(line)
    print(allData)
    
    line1 = allData[0]
    genome = line1[2].lower()

    genomeFolderName = '{}/{}'.format(
        folderName,
        genome
        )

    formatFolder(genomeFolderName,create=True)
    if hub_name == '':
        hub_name = analysis_name

    hubTxt = '.hub.txt'
    hubFileName = '{}/{}{}'.format(
        folderName,
        hub_name,
        hubTxt
        )
    #This is the file you will link to the genome browser that points to the other files in your track hub configuration
    hubFile = open(hubFileName,'w')
    #This is the name of the hub
    hub_line = '{} {}'.format(
        'hub',
        hub_name
        )
    hubFile.write(hub_line + '\n')

    #this is a short label description of this hub
    if hub_short_lab == '':
        hub_short_lab_line = 'shortLabel ' + hub_name
    else:
        hub_short_lab_line = 'shortLabel ' + hub_short_lab
    hubFile.write(hub_short_lab_line+'\n')

    #this is a long label description of this hub
    if hub_long_lab == '':
        hub_long_lab_line = 'longLabel ' + hub_name
    else:
        hub_long_lab_line = 'longLabel ' + hub_long_lab
    hubFile.write(hub_long_lab_line+'\n')

    #points to genomes file
    hubFile.write('genomesFile ' + analysis_name + '.genomes.txt' + '\n')

    #this is the email for the person who set up this track hub
    if EMAIL != '':
        hubFile.write('email ' + EMAIL + '\n')
    else:
        EMAIL = 'email@bcm.edu'
        hubFile.write('email ' + EMAIL + '\n')
    
    hubFile.close()

    genomesFileName = folderName + '/' + analysis_name + '.genomes.txt'
    print(genomesFileName)
    #This is the file that includes all genomes used in this track hub
    genomesFile = open(genomesFileName,'w')

    genomesFile.write('genome ' + genome + '\n')
    genomesFile.write('trackDb ' + genome+'/trackDb.txt' + '\n')

    genomesFile.close()

    #This is the trackDb file which contains information about the files you are going to actually be visualizing

    trackDbFileName = genomeFolderName + '/trackDb.txt'

    trackDbFile = open(trackDbFileName,'w')
    
    if scaled == True:
        wig_str = '_scaled.wig.gz'
    else:
        wig_str = '_treat_afterfiting_all.wig.gz'

    for line in allData:
        name = line[3]
        print('name: ' + name)
        input_wig = wiggle_dir+ name + wig_str
        print('input_wig: ' + input_wig)
        bigwig_name = name + '.bw'
        bigwig_out = '{}/{}'.format(
            genomeFolderName,
            bigwig_name
            )
        if os.path.isfile(input_wig) == True:
            bw_cmd = '{} {} {} {}'.format(
                'wigToBigWig -clip',
                input_wig,
                chrom_sizes,
                bigwig_out
                )
        
            os.system(bw_cmd)

            trackDbFile.write('track ' + name + '\n')
            URL = '{} {}{}{}{}/{}/{}'.format(
                'bigDataUrl',
                urlBase,
                group,
                track_folder,
                analysis_name,
                genome,
                bigwig_name
                )

            trackDbFile.write(URL + '\n')
            trackDbFile.write('shortLabel '+ name + '\n')
            trackDbFile.write('longLabel '+ name + '\n')
            trackDbFile.write('type '+ fileType + '\n')
            trackDbFile.write('color ' + col + '\n\n')

    trackDbFile.close()



#==========================================================================
#===============================GSEA ANALYSIS==============================
#==========================================================================
def wrapGSEA(gctPath,clsPath,sample_1,sample_2,analysis_name,output_folder,metric,gmxPath='',gseaPath='',launch=True):
    '''
    wraps GSEA, creates bash script, and sets up output folder
    Metrics:
        - metric possibilities include Signal2Noise, tTest, Ratio_of_Classes, Diff_of_Classes, log2_Ratio_of_Classes
        - tTest metric requires multiple columns to work
    .gct and .cls files are generated via the normalizeRNAseq.R script in the pipeline, but they can easily be generated.
        GCT files:
          http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT
        CLS files:
          http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CLS
    Gene sets and .gmt files can be generated manually to inspect specific gene sets. The default is a curated gene set of all human genes.
        GMT files:
          http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GMT

    IF YOU GET A STRANGE ERROR, MAKE SURE THE GENES IN YOUR GCT FILE ARE CAPITALIZED.
    '''
    if gmxPath=='':
        #default curated gene set
        gmxPath='/storage/cylin/grail/annotations/gsea/c2.all.v5.1.symbols.gmt'
    if gseaPath=='':
        #default gsea version
        gseaPath = '/storage/cylin/home/cl6/gsea2-3.0_beta_2.jar'

    #create and open bash file  
    gseaBashFilePath = cls_path = '{}{}_gsea.sh'.format(output_folder,analysis_name)
    gseaBashFile = open(gseaBashFilePath,'w')

    #shebang
    gseaBashFile.write('#!/usr/bin/bash\n\n')

    gseaBashFile.write('#COMMAND LINE GSEA CALLS FOR {}\n\n'.format(analysis_name))

    #point where the analysis will live
    gseaOutputFolder = utils.formatFolder('{}{}_GSEA/'.format(output_folder,analysis_name),True)

    #name report
    rptLabel = '{}'.format(analysis_name)

    #generate GSEA command
    gseaCmd_all = 'java -Xmx4000m -cp {} xtools.gsea.Gsea -res {} -cls {}#{}_versus_{} -gmx {} -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label {} -metric {} -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out {} -gui false'.format(gseaPath,gctPath,clsPath,sample_1,sample_2,gmxPath,rptLabel,metric,gseaOutputFolder)

    #write command
    gseaBashFile.write(gseaCmd_all)
    gseaBashFile.close()
    print('writing GSEA output to {}'.format(gseaBashFilePath))
    return gseaBashFilePath
    ##if you want to auto run
    if launch == True:
        os.system('bash {}'.format(gseaBashFilePath))    


