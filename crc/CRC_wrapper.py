
#!/usr/bin/python

'''
The MIT License (MIT)

Copyright (c) 2016 Charles Y. Lin

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

#Description:

#This is the main script for dynamic network comparison
#let's make sure it's python3 compliant

#================================================================================
#=============================DEPENDENCIES=======================================
#================================================================================

import sys
import argparse
import string
import os
from collections import defaultdict


#importing utils package
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = '%s/' % (string.replace(whereAmI,'crc/','')) # need to set this to where this code is stored
sys.path.append(pipeline_dir)


import utils
import pipeline_dfci


whereAmI = os.path.dirname(os.path.realpath(__file__))
#================================================================================
#============================GLOBAL PARAMETERS===================================
#================================================================================

#add locations of files and global parameters in this section


#eventually we will create a crc_wrapper module to hide all of the functions and classes

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

#write your specific functions here

def printAndWriteLog(output,log_file = ''):
    
    '''
    quick helper function that prints and also writes to a log file if the log output is flagged in the options
    '''

    print(output)
    
    if len(log_file) >0:
        try:
            log = open(log_file,'r')
            log.close()
            log = open(log_file,'a')
        except IOError:
            log = open(log_file,'w')

        log.write(output+'\n')
        log.close()



def loadGenome(genome_build):

    '''
    loads annotation for a genome into a genome object
    '''

    #this nested dictionary has all of the useful information and likely will have to be
    #edited so it can be configured any time
    genomeDict = {
        'HG19':{'annot_file':'%sannotation/hg19_refseq.ucsc' % (pipeline_dir),
                'genome_directory':'%sgenomes/' % (pipeline_dir),
                'tf_file':'%scrc/annotation/TFlist_NMid_hg19.txt' % (pipeline_dir),     
                'mask_file':'%sgenomes/Homo_sapiens/UCSC/hg19/Annotation/Masks/hg19_encode_blacklist.bed' % (pipeline_dir),
                }
        }
    
    if genome_build not in genomeDict:
        print('ERROR: UNSUPPORTED GENOME BUILD %s. EXITING NOW' % (genome_build))
        sys.exit()

    #now attempt to load the genome
    genome = Genome(genome_build,genomeDict[genome_build]['genome_directory'],genomeDict[genome_build]['annot_file'])

    #adding additional optional features
    genome.addFeature('tf_file',genomeDict[genome_build]['tf_file'])
    genome.addFeature('mask',genomeDict[genome_build]['mask_file'])


    return genome


def getActivity(genome,dataFile,group1_list,group2_list,output_folder,cutOff = 0.5):
    '''
    gets a list of genes bound by the factor above a cutoff
    '''

    return




def launchMetaRose(group_name,group_list,meta_rose_folder,genome,data_file,stitch,tss):

    '''
    launches meta rose
    '''

    project_folder = utils.formatFolder(os.path.abspath(utils.getParentFolder(data_file)),False)
    macs_folder = '%smacsEnriched/' % (project_folder) #quick hack to make sure input beds are in the right place
    dataDict = pipeline_dfci.loadDataTable(data_file)

    meta_rose_output = utils.formatFolder(meta_rose_folder + group_name,True)

    genome_build = genome.name()
    
    input_string = ','.join(['%s%s' % (macs_folder,dataDict[name]['enrichedMacs']) for name in group_list])
    bam_string = ','.join([dataDict[name]['bam'] for name in group_list])

    meta_cmd = 'python %sROSE2_META.py -g %s -i %s -r %s -o %s -n %s -t %s' % (pipeline_dir,genome_build,input_string,bam_string,meta_rose_output,group_name,tss)
    if stitch != None:
        meta_cmd += ' -s %s' % (stitch)

    #adding a mask if necessary
    if genome.hasFeature('mask'):
        meta_cmd += ' --mask %s' % (genome.returnFeature('mask'))

    bash_path ='%s%s_meta_rose.sh' % (meta_rose_output,group_name)
    bash_file = open(bash_path,'w')
    bash_file.write('#!/usr/bin/bash\n\n')

    bash_file.write('cd %s\n\n' % (pipeline_dir))
    bash_file.write(meta_cmd)
    bash_file.write('\n\n')

    bash_file.close()
    
    print('Wrote META_ROSE command for %s to %s' % (group_name,bash_path))
    print('Launching META_ROSE')
    os.system('bash %s' % (bash_path))





def launchDynamicRose(analysis_name,group1_name,group2_name,group1_list,group2_list,meta_rose_folder,dynamic_rose_folder,genome,data_file,activity_path,inputGFF):

    '''
    launches meta rose
    '''

    project_folder = utils.formatFolder(os.path.abspath(utils.getParentFolder(data_file)),False)

    dataDict = pipeline_dfci.loadDataTable(data_file)

    genome_build = genome.name()
    
    meta_rose_output_1 = utils.formatFolder(meta_rose_folder + group1_name,True)
    meta_rose_output_2 = utils.formatFolder(meta_rose_folder + group2_name,True)

    meta_rose_string = ','.join([meta_rose_output_1,meta_rose_output_2])

    #setting the output
    dynamic_rose_folder = utils.formatFolder(dynamic_rose_folder,True)

    group1_string = ','.join(group1_list)
    group2_string = ','.join(group2_list)
    
    dynamic_cmd = 'python %sdynamicEnhancer_meta.py -g %s -d %s -r %s -o %s --group1 %s --group2 %s --name1 %s --name2 %s -p -m' % (pipeline_dir,genome_build,data_file,meta_rose_string,dynamic_rose_folder,group1_string,group2_string,group1_name,group2_name)
    
    if len(inputGFF) > 0:
        dynamic_cmd += ' --input %s' % (inputGFF)

    bash_path ='%s%s_dynamic_meta.sh' % (dynamic_rose_folder,analysis_name)
    bash_file = open(bash_path,'w')
    bash_file.write('#!/usr/bin/bash\n\n')

    bash_file.write('cd %s\n\n' % (pipeline_dir))
    bash_file.write(dynamic_cmd)
    bash_file.write('\n\n')

    bash_file.close()
    
    print('Wrote DYNAMIC_META command for %s to %s' % (analysis_name,bash_path))
    print('Launching DYNAMIC_META_ROSE')
    os.system('bash %s' % (bash_path))


def launchCRC(data_file,genome,dynamic_rose_output,group_name,group_list,crc_folder,activity_path):
    
    '''
    launches CRC analysis on all bams in a group w/ subpeaks
    #how do we get subpeaks piped through?
    '''
    dataDict = pipeline_dfci.loadDataTable(data_file)    
    bam_string = ','.join([dataDict[name]['bam'] for name in group_list])

    #set up the crc command
    crc_cmd = 'python CRC3.py -e %s -b %s -g %s -o %s -n %s' % (dynamic_rose_output,bam_string,genome.name(),crc_folder,group_name)

    if len(activity_path) > 0:
        crc_cmd += ' --activity %s' % (activity_path)
    bash_path = '%s%s_crc.sh' % (crc_folder,group_name)
    bash_file = open(bash_path,'w')

    bash_file.write('#!/usr/bin/bash\n\n')

    bash_file.write('cd %s\n\n' % (whereAmI))
    bash_file.write(crc_cmd)
    bash_file.write('\n\n')

    bash_file.close()
    
    print('Wrote CRC command for %s to %s' % (group_name,bash_path))
    print('Launching CRC')
    os.system('bash %s' % (bash_path))

    

    



def findCanidateTFs(genome, enhancer_gff, expressedNM, expressionDictNM,
                    bamFile, TFlist, refseqToNameDict, projectFolder, projectName, promoter):
    '''                                                           
    Assign each Super-Enhancer to the closest active TSS to its center
    Return a dictionary keyed by TF that points to a list of loci 
    '''
    
    #loading in the enhancer gff regions
    enhancer_collection = utils.gffToLocusCollection(enhancer_gff)
    enhancer_loci = enhancer_collection.getLoci()


    #loading in the genome and TF info
    annot_file = genome.returnFeature('annot_file')
    startDict = utils.makeStartDict(annot_file)    

    tf_table = utils.parseTable(genome.returnFeature('tf_file'),'\t')
    refID_list = [line[0] for line in tf_table] #creates a list of all NM IDs for TFs

    #make a collection of all TF TSSs
    tssLoci = []
    for refID in refID_list:
        tssLoci.append(utils.makeTSSLocus(refID,startDict,0,0)) #this is a precise 1 coordinate TSS locus
    tssCollection = utils.LocusCollection(tssLoci,50)    



    enhancerTable = [['ENHANCER_ID','CHROM','START','STOP','GENE_LIST']]

    gene_to_enhancer_dict = defaultdict(list)
    # Loop through enhancers
    #all gene nnames stored by refID
    for enhancer in enhancer_loci:
        

        # If the enhancer overlaps a TSS, save it
        overlapping_loci = tssCollection.getOverlap(enhancer, 'both')
        overlapping_refIDs =[locus.ID() for locus in overlapping_loci]

        # Find all gene TSS within 100 kb
        proximal_loci = tssCollection.getOverlap(utils.makeSearchLocus(enhancer,100000,100000),'both')
        proximal_refIDs =[locus.ID() for locus in proximal_loci]
        
        # If no genes are within 100 kb, find the closest active gene within 1 million bp
        closest_refID = []
        if len(overlapping_refIDs) == 0 and len(proximal_refIDs) == 0:
        
            distal_loci = tssCollection.getOverlap(utils.makeSearchLocus(enhancer,1000000,1000000),'both')
            distal_refIDs =[locus.ID() for locus in distal_loci]

            enhancerCenter = (int(enhancer.start()) + int(enhancer.end())) / 2
            distance_list = [abs(enhancerCenter - startDict[geneID]['start'][0])
                             for geneID in distal_refIDs]
            if len(distance_list) > 0:
                closest_refID = [distalGenes[distance_list.index(min(distance_list))]]

        #now we have all potential gene cases
        all_refIDs = overlappingGenes + proximalGenes + closest_refID
        
        #now we get all names and refIDs
        all_refIDs = utils.uniquify([refID for refID in all_refIDs if len(refID) > 0 ])
        all_names = utils.uniquify([startDict[refID]['name'] for refID in all_refIDs])
        
        #first do enhancer level assignment
        names_string = ','.join(all_names)
        enhancer_table.append([enhancer.ID(),enhancer.chr(),enhancer.start(),enhancer.end(),names_string])

        #now do gene level assignment
        for refID in all_refIDs:
            gene_to_enhancer_dict[refID].append(enhancer.ID())

        #an enhancer can be assigned to multiple genes
        #a promoter can only be assigned to 1 gene
        #promoters don't have enhancerIDs so don't add them yet
        #this should just be an enhancer level table
        #followed by a gene level table



        overlappingGenes = utils.uniquify(overlappingGenes)
        proximalGenes = utils.uniquify(proximalGenes)
        for refID in overlappingGenes:
            if proximalGenes.count(refID) == 1:
                proximalGenes.remove(refID)
 

        # If a TSS overlaps an enhancer, assign them together
        if overlappingGenes:
            for gene in overlappingGenes:
                if gene in tf_list:
                    TFtoEnhancerDict[gene].append(enhancer)
                    enhancerAssignment.append([gene, enhancer.chr(), enhancer.start(), enhancer.end(), enhancer.ID()])
                
        # Otherwise, assign the enhancer to the most active gene in 100 kb
        elif not overlappingGenes and proximalGenes:
            highestGene = ''
            highestActivity = 0
            for gene in proximalGenes:
                if expressionDictNM[gene] > highestActivity:
                    highestActivity = expressionDictNM[gene]
                    highestGene = gene
            if highestGene in TFlist:
                TFtoEnhancerDict[gene].append(enhancer)
                enhancerAssignment.append([gene, enhancer.chr(), enhancer.start(), enhancer.end(), enhancer.ID()])
            
        elif not overlappingGenes and not proximalGenes and closestGene:
            if closestGene in TFlist:
                gene = closestGene
                TFtoEnhancerDict[gene].append(enhancer)
                enhancerAssignment.append([gene, enhancer.chr(), enhancer.start(), enhancer.end(), enhancer.ID()])

    # Add promoter is it's not contained in the super
    if promoter:
        for gene in TFtoEnhancerDict.keys():
            promoter = utils.Locus(startDict[gene]['chr'], int(startDict[gene]['start'][0]) - 2000, 
                                   int(startDict[gene]['start'][0]) + 2000, startDict[gene]['sense'])
            overlapBool = False
            for enhancer in TFtoEnhancerDict[gene]:
                if promoter.overlaps(enhancer):
                    overlapBool = True
            if not overlapBool:
                TFtoEnhancerDict[gene].append(promoter)

    seAssignmentFile = projectFolder + projectName + '_ENHANCER_ASSIGNMENT.txt'
    utils.unParseTable(enhancerAssignment, seAssignmentFile, '\t')

    return TFtoEnhancerDict



    
#================================================================================
#===============================MAIN RUN=========================================
#================================================================================


#write the actual script here


def main():

    '''
    using argparse

    '''

    parser = argparse.ArgumentParser(usage='%(prog)s -i DATAFILE -1 GROUP1_NAMES -2 GROUP2_NAMES')

    # required flags
    parser.add_argument("-d", "--data_table", dest="data_table", type=str,
                      help="input a data table with all datasets to be analyzed", required=True)
    parser.add_argument("-1", "--group1", dest="group1", type=str,
                      help="input a comma separated list of all datasets in group1", required=True)
    parser.add_argument("-2", "--group2", dest="group2", type=str,
                      help="input a comma separated list of all datasets in group2", required=True)


    #optional input override
    parser.add_argument("-i", "--input", dest="input", type=str,
                        help="input a gff of regions to analyze", required=False)


    #optional arguments
    parser.add_argument("-n", "--name", dest="name", type=str,
                      help="specify a name for the analysis. Default is drawn from the data table name", required=False)

    parser.add_argument("--group1-name", dest="group1_name", default='GROUP1',type=str,
                      help="Enter a name for group1.  Default is 'GROUP1'", required=False)
    parser.add_argument("--group2-name", dest="group2_name", default='GROUP2',type=str,
                      help="Enter a name for group2.  Default is 'GROUP2'", required=False)

    parser.add_argument("-a", "--activity", dest="activity", type=str,default='',
                      help="a table with active gene names in the first column", required=False)
    parser.add_argument("-t", "--tss", dest="tss", type=int,default=2500,
                      help="Specify a TSS exclusion distance. Default is 2500", required=False)
    parser.add_argument("-s", "--stitch", dest="stitch", type=int,default=None,
                      help="Specify a stitching distance. Default is auto stitching", required=False)



    parser.add_argument("-o", "--output", dest="output", default='./',type=str,
                      help="Enter the output folder. Default is the current working directory", required=False)

    parser.add_argument("--log", dest="log", default='',type=str,
                      help="Enter a path to log output", required=False)



#     # DEBUG OPTION TO SAVE TEMP FILES
#     parser.add_argument("--scale", dest="scale", default='',
#                       help="Enter a comma separated list of scaling factors for your bams. Default is none")
#     parser.add_argument("--save-temp", dest="save", action='store_true', default=False,
#                       help="If flagged will save temporary files made by bamPlot")
#     parser.add_argument("--bed", dest="bed",
#                       help="Add a space-delimited list of bed files to plot")
#     parser.add_argument("--multi-page", dest="multi", action='store_true', default=False,
#                       help="If flagged will create a new pdf for each region")

    args = parser.parse_args()



    #now we can begin to parse the arguments
    
    #=====================================================================================
    #===============================I. PARSING ARGUMENTS==================================
    #=====================================================================================
    #pulling in the data table
    data_file = os.path.abspath(args.data_table)
    dataDict = pipeline_dfci.loadDataTable(data_file)

    #setting naming conventions
    if not args.name:
        analysis_name = data_file.split('/')[-1].split('.')[0]
    else:
        analysis_name = args.name

    #getting the optional input gff
    if args.input:
        inputGFF = args.input
    else:
        inputGFF = ''

    #getting group names
    group1_name = args.group1_name
    group2_name = args.group2_name

    #getting group1 
    group1_string = args.group1
    group1_list = [name for name in string.split(group1_string,',') if len(name) > 0]

    #getting group2
    group2_string = args.group2
    group2_list = [name for name in string.split(group2_string,',') if len(name) > 0]

    #checking that all datasets are in the data table
    for name in group1_list + group2_list:
        if name not in dataDict:
            print('ERROR: DATASET %s NOT FOUND IN DATA TABLE %s. EXITING NOW' % (name,data_file))
            sys.exit()

    #loading in the genome object from the data table
    genome_list = utils.uniquify([dataDict[name]['genome'] for name in group1_list + group2_list])
    if len(genome_list) > 1:
        print('ERROR: ATTEMPTING TO ANALYZE DATASETS FROM MULTIPLE GENOME BUILDS. EXITING NOW.')
        sys.exit()

    
    #the load genome function has an assertion test to make sure the genome is supported
    genome = loadGenome(genome_list[0])

    
    parent_folder = utils.formatFolder(args.output,True)
    output_folder = utils.formatFolder(parent_folder + analysis_name,True)


    #these are the user defined optional arguments
    tss = int(args.tss)

    stitch = args.stitch
    print('stitch')
    print(stitch)

    
    #list of active genes to constrain analysis 
    if len(args.activity) == 0:
        #assumes all genes are active unless told otherwise
        #activity_path,activity_table = getActivity() # fix this function
        print('using all active genes')
    else:
        activity_path = args.activity
        activity_table = utils.parseTable(activity_path,'\t')




    print('\n\n#======================================\n#===========I. DATA SUMMARY============\n#======================================\n')

    print('Analyzing datasets described in %s\n' % (data_file))

    print('Name for the analysis: %s\n' % (analysis_name))
    print('Using genome: %s\n' % (genome.name()))


    
    print('%s datasets: %s\n' % (group1_name,group1_string))
    print('%s datasets: %s\n' % (group2_name,group2_string))

    if len(activity_path) > 0:
        print('Identified %s active genes in the analysis using %s as a list of active genes' % (len(activity_table),activity_path))
    else:
        print('Identified %s active genes in the analysis using aggregate data from %s and %s' % (len(activity_table),group1_name,group2_name))
    print('Writing output to: %s\n' % (output_folder))


    #=====================================================================================
    #======================II. DEFINING CIS-REGULATORY ELEMENTS===========================
    #=====================================================================================


    print('\n\n#======================================\n#=II. MAPPING CIS-REGULATORY ELEMENTS==\n#======================================\n')



    #crc_wrapper will act at the group level and not consider individual datasets
    #since a data table is used as the input, the code will rely heavily on pipeline_dfci
    #embedded tools

    #1. first we need to run meta rose using default parameters and check the output
    #exists for each group

    meta_rose_folder = utils.formatFolder(output_folder + 'meta_rose/',True)

    group1_output = '%s%s/%s_AllEnhancers.table.txt' % (meta_rose_folder,group1_name,group1_name)

    group2_output = '%s%s/%s_AllEnhancers.table.txt' % (meta_rose_folder,group2_name,group2_name)
    #print(group1_output)
    #print(group2_output)

    #for each output check to see if they exist
    #if not launch

    try:
        foo = open(group1_output,'r')
    except IOError:
        print('No META_ROSE output found for %s. Running META_ROSE now' % (group1_name))
        launchMetaRose(group1_name,group1_list,meta_rose_folder,genome,data_file,stitch,tss)
        
    try:
        foo = open(group2_output,'r')
    except IOError:
        print('No META_ROSE output found for %s. Running META_ROSE now' % (group2_name))
        launchMetaRose(group2_name,group2_list,meta_rose_folder,genome,data_file,stitch,tss)



    #now check for completion
    if utils.checkOutput(group1_output,1,10):
        print('META_ROSE finished for %s' % (group1_name))
    else:
        print('META_ROSE timed out for %s. EXITING NOW.' % (group1_name))
        sys.exit()

    if utils.checkOutput(group2_output,1,10):
        print('META_ROSE finished for %s' % (group2_name))
    else:
        print('META_ROSE timed out for %s. EXITING NOW.' % (group2_name))
        sys.exit()


    #Meta rose does not give all regions that are SE in at least one sample
    #and can be blown out by amplicons etc...
    #sooo we need to run clustering to generate a good input gff
    #ideally we just rewrite dynamic meta to run off of clustering output
    #until we do that let's just overwrite w/ an input gff
    

    print('Comparing cis-regulatory landscapes of %s and %s' % (group1_name,group2_name))
    dynamic_rose_folder = utils.formatFolder(output_folder + 'dynamic_meta_rose/',True)

    #here we will use the rank table as the primary output
    dynamic_rose_output = '%soutput/%s_%s_%s_merged_MERGED_SUPERS_RANK_TABLE.txt' % (dynamic_rose_folder,genome.name(),group1_name,group2_name)
    
    try:
        foo = open(dynamic_rose_output,'r')
    except IOError:
        print('No DYNAMIC_ROSE output found for %s. Running DYNAMIC_ROSE now' % (analysis_name))
        launchDynamicRose(analysis_name,group1_name,group2_name,group1_list,group2_list,meta_rose_folder,dynamic_rose_folder,genome,data_file,activity_path,inputGFF)

    if utils.checkOutput(dynamic_rose_output,1,10):
        print('DYNAMIC_ROSE finsihed for %s' % (analysis_name))
    else:
        print('DYNAMIC_ROSE analysis timed out for %s. EXITING NOW.' % (analysis_name))
        sys.exit()




    #=====================================================================================
    #======================III. IDENTIFYING TF NODES IN NETWORK===========================
    #=====================================================================================


    print('\n\n#======================================\n#===III. RUNNING CIRCUITRY ANALYSIS====\n#======================================\n')




    #now we want to call circuitry on each group... ok to have different subpeaks and motif calls
    #if as a first approximation we weight by the overall enhancer




    crc_folder = utils.formatFolder('%scrc/' % (output_folder),True)



    #for all
    all_crc_folder = utils.formatFolder('%s%s' % (crc_folder,analysis_name),True)
    launchCRC(data_file,genome,dynamic_rose_output,analysis_name,group1_list+group2_list,all_crc_folder,activity_path)



    #for group1
    group1_crc_folder = utils.formatFolder('%s%s' % (crc_folder,group1_name),True)
    launchCRC(data_file,genome,dynamic_rose_output,group1_name,group1_list,group1_crc_folder,activity_path)

    #for group2
    group2_crc_folder = utils.formatFolder('%s%s' % (crc_folder,group2_name),True)
    launchCRC(data_file,genome,dynamic_rose_output,group2_name,group2_list,group2_crc_folder,activity_path)







if __name__=="__main__":
    main()
