#!/usr/bin/python

import sys
import os
import argparse
import string
import utils

def split_sam(suffix,samfile,header1,header2,outfile1,outfile2):

    '''
    splits the sam by suffix
    cats on the header
    '''
    print(suffix)
    temp_outfile1 = outfile1+'.temp'
    temp_outfile2 = outfile2+'.temp'
    
    temp1 = open(temp_outfile1,'w')
    temp2 = open(temp_outfile2,'w')

    sam = open(samfile,'r')
    ticker =0
    print('lines processed')
    for line in sam:

        ticker+=1
        if ticker%1000000 == 0:
            print(ticker)
        if line[0] == '@':
            continue

        if line.split('\t')[2][-4:] == suffix:
            line = string.replace(line,suffix,'')
            temp2.write(line)
        else:
            temp1.write(line)

    temp1.close()
    temp2.close()

    cmd1 = 'cat %s %s > %s' % (header1,temp_outfile1,outfile1)
    cmd2 = 'rm %s' % (temp_outfile1)


    cmd3 = 'cat %s %s > %s' % (header2,temp_outfile2,outfile2)
    cmd4 = 'rm %s' % (temp_outfile2)
    
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)

    #now convert to bam
    print('CONVERTING TO BAM')
    bamfile1 = string.replace(outfile1,'.sam','.bam')
    bamfile2 = string.replace(outfile2,'.sam','.bam')
    cmd5 = 'samtools view -bS %s > %s' % (outfile1,bamfile1)
    cmd6 = 'samtools view -bS %s > %s' % (outfile2,bamfile2)
    
    sorted_bamfile1_prefix = string.replace(bamfile1,'.bam','.sorted')
    sorted_bamfile1 = sorted_bamfile1_prefix + '.bam'

    sorted_bamfile2_prefix = string.replace(bamfile2,'.bam','.sorted')
    sorted_bamfile2 = sorted_bamfile2_prefix + '.bam'


    cmd7 = 'samtools sort %s %s' % (bamfile1,sorted_bamfile1_prefix)
    cmd8 = 'samtools sort %s %s' % (bamfile2,sorted_bamfile2_prefix)

    cmd9 = 'samtools index %s' % (sorted_bamfile1)
    cmd10 = 'samtools index %s' % (sorted_bamfile2)

    os.system(cmd5)
    os.system(cmd6)

    print('SORTING AND INDEXING BAMS')
    os.system(cmd7)
    os.system(cmd8)
    os.system(cmd9)
    os.system(cmd10)

    print('DELETING TEMP FILES')

    cmd11 = 'rm %s' % (temp_outfile1)
    cmd12 = 'rm %s' % (temp_outfile2)
    cmd13 = 'rm %s' % (bamfile1)
    cmd14 = 'rm %s' % (bamfile2)

    os.system(cmd11)
    os.system(cmd12)
    os.system(cmd13)
    os.system(cmd14)




def main():
    """
    main run function
    """

    #usage = "usage: %prog [options] -g [GENOME] -b [SORTED BAMFILE(S)] -i [INPUTFILE] -o [OUTPUTFOLDER]"
    parser = argparse.ArgumentParser(usage='%(prog)s -i SAMFILE -g REF_GENOME -s SPIKE_GENOME')

    # required flags
    parser.add_argument("-i", "--input", dest="inputSam", type=str,
                        help="Enter a sam file", required=False)
    parser.add_argument("-g", "--genome", dest="genome", type=str,
                        help="specify the main reference genome", required=False)
    parser.add_argument("-s", "--spike", dest="spike", type=str,
                        help="specify the spike in  genome", required=False)


    parser.add_argument("-d", "--dest", dest="dest", type=str,
                        help="specify an optional destination for the final bams to move to", required=False)

    args = parser.parse_args()

    print(args)
    
    if args.inputSam and args.genome and args.spike:

        print('FORMATTING %s FOR CHIP_RX USING REFERENCE GENOME %s and SPIKE_IN GENOME %s' % (args.inputSam,args.genome,args.spike))
        samPath = args.inputSam

        if string.upper(samPath).count('.SAM') == 0:
            print('ERROR, file must end in .sam or .SAM')
            sys.exit()

        
        #get the headers
        genome_string = string.upper('%s_%s' % (args.genome,args.spike))

        
        genomeDict = {'RN6_DM6':['/grail/genomes/Rattus_norvegicus/UCSC/rn6/Sequence/Bowtie2Index_dm6/rn6_header.SAM','/grail/genomes/Rattus_norvegicus/UCSC/rn6/Sequence/Bowtie2Index_dm6/dm6_header.SAM'],
                      }

        if genomeDict.has_key(genome_string) == False:
            print('ERROR: UNSUPPORTED GENOME COMBINATION %s' % (genome_string))
            sys.exit()
        else:
            print('IDENTIFIED HEADERS FOR GENOME COMBINATION %s' %(genome_string))


        header1 = genomeDict[genome_string][0]
        header2 = genomeDict[genome_string][1]
        suffix = string.lower('_%s' % (args.spike))

        outfile1 = string.replace(samPath,samPath[-4:],'.%s%s' % (args.genome,samPath[-4:]))
        outfile2 = string.replace(samPath,samPath[-4:],'.%s%s' % (args.spike,samPath[-4:]))
        split_sam(suffix,samPath,header1,header2,outfile1,outfile2)

        #move stuff to destination folder
        if args.dest:
            bamFolder = utils.formatFolder(args.dest,False)

            samFolder = utils.getParentFolder(samPath)

            mv_cmd = 'mv %s*bam* %s' % (samFolder,bamFolder)
            print('MOVING BAMS FROM %s TO %s' % (samFolder,bamFolder))
            os.system(mv_cmd)


    else:
        parser.print_help()
        sys.exit()

if __name__ == "__main__":
    main()


