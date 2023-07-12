#!/usr/bin/python


'''
The MIT License (MIT)

Copyright (c) 2018 Charles Lin

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


#Functions for calculating delta in and out degree from CRC output edge tables
#First commit with some hard paths


#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================

#requires bamliquidator and the linlab pipeline installed
#runs on an edge table from CRC output

import sys, os
# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = whereAmI.replace('/crc','')
print(pipeline_dir)

sys.path.append(whereAmI)
sys.path.append(pipeline_dir)

import pipeline_dfci
import utils
import string
import numpy
from scipy import stats
import os
import re
from collections import defaultdict
#==========================================================================
#============================PARAMETERS====================================
#==========================================================================

#hard coded parameters can go here during debug

#crc_folder <- standard CRC output
#chip_dataFile <- a linlab format data table w/ chip data to be mapped at edges
#analysis_name <- analysis name used in CRC. Typical CRC output files will include [ANALYSIS_NAME]_EDGE_TABLE.txt e.g. if NIBR_EDGE_TABLE.txt then analysis_name = NIBR
#group1_list,group2_list <- names of datasets in each group (specified in the chip_dataFile)
#output = path to write to. default will write to the CRC folder

#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================
#Example data file for a chip dataset
#ChIP-Seq
chip_dataFile = '%sdata_tables/NIBR_CHIP_TABLE.txt' % (projectFolder)





        
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~CALCULATING CHANGES IN BRD4 OUT DEGREE BY TF EDGES~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def tf_edge_delta_out(crc_folder,chip_dataFile,analysis_name,group1_list,group2_list,output=''):

    '''
    calculates changes in brd4 out degree at each predicted motif occurrence this is by subpeaks
    '''
    crc_folder = utils.formatFolder(crc_folder,False)
    edge_path = '%s%s_EDGE_TABLE.txt' % (crc_folder,analysis_name)

    #make a gff of the edge table
    edge_table = utils.parseTable(edge_path,'\t')
    edge_gff = []
    for line in edge_table[1:]:
        gff_line = [line[2],'%s_%s' % (line[0],line[1]),'',line[3],line[4],'','.','','%s_%s' % (line[0],line[1])]
        edge_gff.append(gff_line)

    edge_gff_path = '%s%s_EDGE_TABLE.gff' % (crc_folder,analysis_name)
    utils.unParseTable(edge_gff,edge_gff_path,'\t')
    
    #direct the output to the crc folder
    signal_path = '%s%s_EDGE_TABLE_signal.txt' % (crc_folder,analysis_name)
    
    #get a list of all chip datasets
    all_chip_list = group1_list + group2_list
    if utils.checkOutput(signal_path,0,0) == False:
        signal_table_list = pipeline_dfci.map_regions(chip_dataFile,[edge_gff_path],mappedFolder,signalFolder,all_chip_list,True,signal_path,extendReadsTo=100)
        print(signal_table_list)
    else:
        print('Found previous signal table at %s' % (signal_path))
    
    #now bring in the signal table as a dictionary using the locus line as the id
    print('making log2 group1 over group2 signal table at edges')
    signal_table = utils.parseTable(signal_path,'\t')
    signal_dict = defaultdict(float)
    
    #figure out columns for group1 and group2
    group2_columns = [signal_table[0].index(name) for name in group2_list]
    group1_columns = [signal_table[0].index(name) for name in group1_list]
    group2_signal_vector = []
    group1_signal_vector = []
    for line in signal_table[1:]:
        group2_signal = numpy.mean([float(line[col]) for col in group2_columns])
        group1_signal = numpy.mean([float(line[col]) for col in group1_columns])
        
        group2_signal_vector.append(group2_signal)
        group1_signal_vector.append(group1_signal)

    group2_median = numpy.median(group2_signal_vector)
    group1_median = numpy.median(group1_signal_vector)

    print('group2 median signal (rpm/bp)')
    print(group2_median)
    print('group1 median signal (rpm/bp)')
    print(group1_median)

    #now that we have the median, we can take edges where at least 1 edge is above the median
    #and both are above zero and generate a new table w/ the fold change 

    signal_filtered_path = string.replace(signal_path,'.txt','_filtered.txt')
    if utils.checkOutput(signal_filtered_path,0,0):
        print('Found filtered signal table for edges at %s' % (signal_filtered_path))
        signal_table_filtered = utils.parseTable(signal_filtered_path,'\t')
    else:
    
        signal_table_filtered = [signal_table[0]+['GROUP2_MEAN','GROUP1_MEAN','LOG2_GROUP1_OVER_GROUP2']]
        for line in signal_table[1:]:
            group2_signal = numpy.mean([float(line[col]) for col in group2_columns])
            group1_signal = numpy.mean([float(line[col]) for col in group1_columns])

            if (group2_signal > group2_median or group1_signal > group1_median) and min(group2_signal,group1_signal) >0:
                delta = numpy.log2(group1_signal/group2_signal)
                new_line = line + [group2_signal,group1_signal,delta]
                signal_table_filtered.append(new_line)

        utils.unParseTable(signal_table_filtered,signal_filtered_path,'\t')

    #now get a list of all TFs in the system
    tf_list = utils.uniquify([line[0].split('_')[0] for line in signal_table_filtered[1:]])
    tf_list.sort()
    print(tf_list)

    out_degree_table = [['TF_NAME','EDGE_COUNT','DELTA_MEAN','DELTA_MEDIAN','DELTA_STD','DELTA_SEM']]
    
    for tf_name in tf_list:
        print(tf_name)
        edge_vector = [float(line[-1]) for line in signal_table_filtered[1:] if line[0].split('_')[0] == tf_name]

        edge_count = len(edge_vector)
        delta_mean = round(numpy.mean(edge_vector),4)
        delta_median = round(numpy.median(edge_vector),4)
        delta_std = round(numpy.std(edge_vector),4)
        delta_sem = round(stats.sem(edge_vector),4)
        tf_out_line = [tf_name,edge_count,delta_mean,delta_median,delta_std,delta_sem]
        out_degree_table.append(tf_out_line)
        
    if output == '':
        #set final output
        output_path = '%s%s_EDGE_DELTA_OUT.txt' % (crc_folder,analysis_name)

    else:
        output_path = output

    utils.unParseTable(out_degree_table,output_path,'\t')
    print(output_path)
    return(output_path)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~CALCULATING CHANGES IN BRD4 IN DEGREE AT TFS~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def tf_edge_delta_in(crc_folder,chip_dataFile,analysis_name,group1_list,group2_list,output=''):

    '''
    calculates changes in BRD4 at in degree edges
    '''
    crc_folder = utils.formatFolder(crc_folder,False)
    enhancer_tf_path = '%s%s_ENHANCER_TF_TABLE.txt' % (crc_folder,analysis_name)

    #make a gff of the edge table
    enhancer_tf_table = utils.parseTable(enhancer_tf_path,'\t')
    enhancer_tf_gff = []
    for line in enhancer_tf_table[1:]:
        gff_line = [line[1],line[0],'',line[2],line[3],'','.','',line[0]]
        enhancer_tf_gff.append(gff_line)

    enhancer_tf_gff_path = '%s%s_ENHANCER_TF_TABLE.gff' % (crc_folder,analysis_name)
    utils.unParseTable(enhancer_tf_gff,enhancer_tf_gff_path,'\t')
    
    #direct the output to the crc folder
    signal_path = '%s%s_ENHANCER_TF_TABLE_signal.txt' % (crc_folder,analysis_name)
    
    all_chip_list = group1_list + group2_list
    if utils.checkOutput(signal_path,0,0) == False:
        signal_table_list = pipeline_dfci.map_regions(chip_dataFile,[enhancer_tf_gff_path],mappedFolder,signalFolder,all_chip_list,True,signal_path,extendReadsTo=100)
        print(signal_table_list)
    else:
        print('Found previous signal table at %s' % (signal_path))


    #now bring in the signal table as a dictionary using the locus line as the id
    print('making an enhancer signal dict')
    signal_table = utils.parseTable(signal_path,'\t')
    group1_signal_dict = defaultdict(float)
    group2_signal_dict = defaultdict(float)
    
    #signal here is calculated as AUC
    
    #figure out columns for group2 and group1
    group2_columns = [signal_table[0].index(name) for name in group2_list]
    group1_columns = [signal_table[0].index(name) for name in group1_list]

    for line in signal_table[1:]:
        region_coords = [int(x) for x in line[1].split(':')[-1].split('-')]
        region_length = region_coords[1] - region_coords[0]
        group2_signal = region_length*numpy.mean([float(line[col]) for col in group2_columns])
        group1_signal = region_length*numpy.mean([float(line[col]) for col in group1_columns])
        group1_signal_dict[line[0]] = group1_signal
        group2_signal_dict[line[0]] = group2_signal

    #now grab the gene table
    gene_tf_path = '%s%s_GENE_TF_TABLE.txt' % (crc_folder,analysis_name)
    gene_tf_table = utils.parseTable(gene_tf_path,'\t')
    group1_tf_dict = defaultdict(float)
    group2_tf_dict = defaultdict(float)
    
    for line in gene_tf_table[1:]:
        group1_tf_dict[line[0]] += group1_signal_dict[line[-1]]
        group2_tf_dict[line[0]] += group2_signal_dict[line[-1]]

    tf_list = utils.uniquify([line[0] for line in gene_tf_table[1:]])
    
    tf_list.sort()
    print(tf_list)
    
    in_degree_table = [['TF_NAME','GROUP1_IN','GROUP2_IN','LOG2_GROUP1_vs_GROUP2']]
    
    for tf_name in tf_list:
        group1_signal = round(group1_tf_dict[tf_name],4)
        group2_signal = round(group2_tf_dict[tf_name],4)
        delta = round(numpy.log2(group1_signal/group2_signal),4)
        new_line = [tf_name,group1_signal,group2_signal,delta]
        in_degree_table.append(new_line)

    
    if output == '':
        #set final output
        output_path = '%s%s_TF_DELTA_IN.txt' % (crc_folder,analysis_name)

    else:
        output_path = output

    utils.unParseTable(in_degree_table,output_path,'\t')
    print(output_path)
    return(output_path)
                  
                

