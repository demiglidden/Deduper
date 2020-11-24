#!/usr/bin/env python

####################################################################################
# This program takes a SAM file sorted by chromosomes and removes PCR duplicates   #
# returning a deduplicated file, a duplicates file, as well as files for           #
# misindexed and unmapped reads. Finally, it will provide a summary file with read #
# counts for each category.                                                        #
####################################################################################


import re
import argparse

def get_args():
    parser = argparse.ArgumentParser("a program to remove PCR duplicates from a SAM file")
    parser.add_argument("-sf", "--samfile", type=str, help="SAM file to pass through program, REQUIRED", required=True)
    parser.add_argument("-u", "--umi", type=str, help="file with known umis, REQUIRED", required=True)
    parser.add_argument("-ns", "--namingscheme", type=str, help="naming scheme to be added at the end of all output file names, REQUIRED", required=True)
    parser.add_argument("-pe", "--pairedend", action="store_true", help="paired-end data not supported", required=False)
    parser.add_argument("-r", "--randomers", action="store_true", help="randomers not supported, script requires known umis", required=False)


    return parser.parse_args()

args = get_args()

###assign arguments to variables for use inside of program
sam = args.samfile
umis = args.umi
namingscheme = args.namingscheme
pe = args.pairedend
rand = args.randomers

if pe:
    raise Exception("This program does not support paired-end data.")

if rand == True:
    raise Exception("This program does not support randomers. Known umis are required. Please pass umis to program using -u or --umi.")

def get_mapped(record):
    '''This function will check to see if the read is mapped or not, field 3 of the SAM file.'''
    record = record.strip()
    x = re.search('\s([0-9]+)(\s.+\s[0-9a-zA-Z]+\s)', record)
    flag = int(x[1])
    if((flag & 4) != 4):
        mapped = True
    else: 
        mapped = False    
    return mapped


def get_strand(record):
    '''This function will check to see if the read is forward or reverse, field 3 of SAM file.'''
    record = record.strip()
    x = re.search('\s([0-9]+)(\s.+\s[0-9a-zA-Z]+\s)', record)
    flag = int(x[1])
    if((flag & 16) != 16):
        strand = True
    else: 
        strand = False       
    return strand


def get_adjstart_fwd(record):
    '''This function will subtract the leftmost softclipping from the alignment start position to 
    get the "true" alignment start position, using the CIGAR string, feild 6 of the SAM file.'''
    record = record.strip()
    x = re.search('(:[A-Z]{8}\s)(\d+\s)(\S+)(\s\d+)(\s\S+)(\S+\s)(\S+\s)', record)
    cigar = str(x[7])
    cigar = cigar.strip()
    y = re.search('(:[A-Z]{8}\s)(\d+\s)(\S+)(\s\d+)', record)
    startpos = int(y[4])
    if 'S' in cigar:
        try:
            clip = int(cigar.split('S')[0])
            start_pos_fwd = startpos - clip
        except:
            start_pos_fwd = startpos
    else:
        start_pos_fwd = startpos
    
    return start_pos_fwd


def get_adjstart_rev(record):
    '''This function will subtract the leftmost softclipping from the alignment start position to 
    get the "true" alignment start position, using the CIGAR string, feild 6 of the SAM file.
    It also accounts for intron gaps (N) and deletions from the reference (D)'''
    record = record.strip()
    x = re.search('(:[A-Z]{8}\s)(\d+\s)(\S+)(\s\d+)(\s\S+)(\S+\s)(\S+\s)', record)
    cigar = str(x[7])
    cigar = cigar.strip()
    m_pos = cigar.index("M")
    cigar_atoms = []
    for x in re.findall('[0-9]+[A-Z]{1}', cigar[m_pos:]):
        if x[-1] == 'M':
            cigar_atoms.append(x)
        elif x[-1] == 'S':
            cigar_atoms.append(x)
        elif x[-1] == 'D':
            cigar_atoms.append(x)
        elif x[-1] == 'N':
            cigar_atoms.append(x)
        else:
            continue
    y = re.search('(:[A-Z]{8}\s)(\d+\s)(\S+)(\s\d+)', record)
    startpos = int(y[4])
    for atom in cigar_atoms:
        numb = re.search('[0-9]+', atom)
        numb = int(numb[0])
        startpos += numb
    start_pos_rev = startpos
    
    return start_pos_rev


def get_umi(record):
    '''This function will extract the UMI from the header, field 1 of the SAM file.'''
    record = record.strip()
    umi = re.search('(:)([A-Z]{8}\s)', record)
    umi = umi[2]
    umi = umi.strip()
    
    return umi


def get_chrom(record):
    '''This function will extract the chromosome from the read, field 3 of the SAM file.'''
    record = record.strip()
    x = re.search('(:[A-Z]{8}\s)(\d+\s)(\S+)', record)
    chrom = x[3]

    return chrom

####ALGORITHM

###open files to write to:

deduped_out = open("deduplicated_"+namingscheme+".sam", "w")
duplicates_out = open("duplicates_"+namingscheme+".sam", "w")
unmapped_out = open("unmapped_"+namingscheme+".sam", "w")
misindexed_out = open("misindexed_"+namingscheme+".sam", "w")
summary_out = open("summary_"+namingscheme+".txt", "w")

###initialize chromosome dictionary
chrom_dict = {}

###read umi file
umis_dict = {}
with open (umis, "r") as uh:
    for line in uh:
        line = line.strip()
        umis_dict[line] = 0

###initialize counters
deduped_count = 0
duplicates_count = 0
unmapped_count = 0
misindexed_count = 0

chrommarker = ""

###read samfile
with open (sam, "r") as fh:
    for record in fh:
        if record[0] == '@':
            deduped_out.write(record)
            duplicates_out.write(record)
            unmapped_out.write(record)
        else:
            mapped = get_mapped(record)
            if mapped == False:
                unmapped_out.write(record)
                unmapped_count += 1
            else:
                chrom = get_chrom(record)
                umi = get_umi(record)
                if chrommarker != chrom:
                    chrommarker = chrom
                    chrom_dict.clear()
                if umi not in umis_dict:
                    misindexed_out.write(record)
                    misindexed_count += 1
                else:
                    strand = get_strand(record)
                    if strand == True:
                        startpos = get_adjstart_fwd(record)
                    else:
                        startpos = get_adjstart_rev(record)
                    #print(startpos)
                    uniqid = str(chrom)+str(umi)+str(strand)+str(startpos)
                    if uniqid not in chrom_dict:
                        chrom_dict[uniqid] = 1
                        deduped_out.write(record)
                        deduped_count += 1
                    else:
                        chrom_dict[uniqid] += 1
                        duplicates_out.write(record)
                        duplicates_count += 1


                

###write summary file
summary_out.write("Deduplicated records:"+str(deduped_count)+"\n"+"Duplicate records:"+str(duplicates_count)+"\n"+"Unmapped records:"+str(unmapped_count)+"\n"+"Misindexed records:"+str(misindexed_count)+"\n")

###closing files
deduped_out.close()
duplicates_out.close()
misindexed_out.close()
unmapped_out.close()
summary_out.close()
