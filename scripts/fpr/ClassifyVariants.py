#!/usr/bin/env python
import argparse
import os, sys
import re
import pandas

pandas.set_option('mode.chained_assignment', None)

parser = argparse.ArgumentParser(description=' ')
parser.add_argument('--bedfile')
parser.add_argument('--countsfile')

args = parser.parse_args()


# file containing the CountReads output for the given cell
countsfile = args.countsfile
bedfile = args.bedfile

# check that expanded_hetsnps bedfile is not empty, if empty, exit bc we have no output
if os.stat(bedfile).st_size == 0:
    print('')
    sys.exit()



counts = open(countsfile,'r')
lines = counts.readlines()

# create table from cell=specific bed file that stores position of each variant for classifying
# bed=pandas.read_table("/gpfs/commons/groups/landau_lab/skao/pta/invitro_comparison/expanded_hetSNPs/scan2phycall.disagreement.C2.expanded.hetSNPs.bed",sep='\t',header=None)
bed=pandas.read_table(bedfile,sep='\t',header=None)
positions = pandas.DataFrame(list(zip(bed.iloc[:,0],bed.iloc[:,1],bed.iloc[:,8],bed.iloc[:,4])), columns=['chr', 'variant', 'hetSNP', 'branch'])
positions['category'] = "NA"
positions['variant'] += 100 # we expanded the bed file intervals by 100 so +100 to get original pos
positions['REF_ALT'] = 0 
positions['REF_REF'] = 0 
positions['ALT_REF'] = 0 
positions['ALT_ALT'] = 0 
positions['num_refalt'] = 0 # number of reads containing a combo of ref/alt alleles 
positions['num_nonmatching'] = 0 # number of nonmatching reads

# create table to store output
output_table = pandas.DataFrame(columns=['chr','variant','hetSNP','branch', 'category', 
                                         'REF_ALT', 'REF_REF', 'ALT_REF', 'ALT_ALT', 
                                         'num_refalt', 'num_nonmatching'])


# iterate over each line in the table containing our variant and hetSNP positions
# and categorize as FP or TP (fill in 'category' column)
truepos = 0
falsepos = 0
for r in positions.index:
    chr = positions['chr'][r]
    varposition = positions['variant'][r]
    hetsnpposition = positions['hetSNP'][r]
    
    # the dataframe to store the number of each type of matching read
    refalt_match_counts = pandas.DataFrame(columns=['REF_ALT', 'REF_REF', 'ALT_REF', 'ALT_ALT'], index=range(1))
    for col in refalt_match_counts.columns:
        refalt_match_counts[col].values[:] = 0
    
    l = 0 # index of current line
    
    #print("Variant ", str(varposition))
    
    for line in lines:
        readcatcount = 0 # read category count; the number of categories of reads found for this var/hetsnp pair
        
        refalt_reads = 0
        nonmatching_reads = 0
        # if line contains varposition and hetsnpposition and next line contains "Variant"
        if re.search(str(varposition), line) and re.search(str(hetsnpposition), line) and re.search("Variant", lines[l + 1]):
            
            # iterate over following lines 
            # and count number of lines that show support for ref/alt alleles 
            i = l + 2 # current line to count read categories
            while i < len(lines) and lines[i].startswith("\tReads supporting"): # only count reads with some combo of alt/ref
                #print(lines[i])
                # update category of count type 
                linestr = lines[i].split() 
                
                if re.search('REF.*ALT', lines[i]):
                    # debugging:
                    #print("REF_ALT added: ", linestr[-1])
                    refalt_match_counts['REF_ALT'][0] += int(linestr[-1])
                if re.search('REF.*REF', lines[i]):
                    # debugging:
                    #print("REF_REF added: ", linestr[-1])
                    refalt_match_counts['REF_REF'][0] += int(linestr[-1])
                if re.search('ALT.*REF', lines[i]):
                    # debugging:
                    #print("ALT_REF added: ", linestr[-1])
                    refalt_match_counts['ALT_REF'][0] += int(linestr[-1])
                if re.search('ALT.*ALT', lines[i]):
                    # debugging:
                    #print("ALT_ALT added: ", linestr[-1])
                    refalt_match_counts['ALT_ALT'][0] += int(linestr[-1])
                
                # add number of refalt reads to our total count
                refalt_reads += int(linestr[-1])
                
                
                readcatcount += 1
                i += 1
            
            # add read type counts to positions
            positions.at[r, 'REF_ALT'] = refalt_match_counts['REF_ALT'][0]
            positions.at[r, 'REF_REF'] = refalt_match_counts['REF_REF'][0]
            positions.at[r, 'ALT_REF'] = refalt_match_counts['ALT_REF'][0]
            positions.at[r, 'ALT_ALT'] = refalt_match_counts['ALT_ALT'][0]
                
            # iterate over lines again and count number of nonmatching reads
            i += 1 # start after the line we ended at above
            while i < len(lines) and lines[i].startswith("\tReads with"): 
                # add number of nonmatching reads to our count
                linestr = lines[i].split()
                nonmatching_reads += int(linestr[-1])
                
                i += 1

                
            
            # fill out table according to how many read categories we observe
            if readcatcount > 2:
                # positions['category'][r] = 'FP'
                positions.at[r, 'category'] = 'FP'
                positions.at[r, 'num_refalt'] = refalt_reads
                positions.at[r, 'num_nonmatching'] = nonmatching_reads
                falsepos += 1
            elif readcatcount <= 2: 
                # positions['category'][r] = 'TP'
                positions.at[r, 'category'] = 'TP'
                positions.at[r, 'num_refalt'] = refalt_reads
                positions.at[r, 'num_nonmatching'] = nonmatching_reads
                truepos += 1
        l += 1
    

# drop duplicates from the table (some variants have more than one nearby hetSNP, we want to only use the closest one)
# for r in positions.index:
#     duplicates = positions[positions.variant == positions['variant'][r]]
#     duplicates['distance'] = abs(positions['hetSNP'] - positions['variant'])
    
#     # if we find duplicates, find the one where the hetSNP is closest to the variant 
#     if len(duplicates) > 1:
#         duplicates = duplicates.loc[duplicates['distance'] == duplicates['distance'].min()]
    
#     # keep only the variant/hetsnp pair with the closest hetsnp
#     # output_table = output_table.append(duplicates)
#     # if the row already exists in the table we cannot add it again
#     output_table = pandas.concat([output_table, duplicates], ignore_index = True, sort=False)
    
# output_table = output_table.drop_duplicates()        
# output_table['variant'] += 1 # add 1 to get original pos


# no longer dropping duplicates
output_table = positions
positions
print(countsfile)
# print variant and its category
truepos_dedup = 0
falsepos_dedup = 0
print('chr:variant\thetSNP\tbranch\tcategory\tREF_ALT\tREF_REF\tALT_REF\tALT_ALT\tmatching_reads\tnonmatching_reads\ttotal_reads')
for r in output_table.index:
    print('%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % 
          (str(output_table['chr'][r]), 
           str(output_table['variant'][r]), 
           str(output_table['hetSNP'][r]),
           str(output_table['branch'][r]), 
           str(output_table['category'][r]),
           str(output_table['REF_ALT'][r]),
           str(output_table['REF_REF'][r]),
           str(output_table['ALT_REF'][r]),
           str(output_table['ALT_ALT'][r]),
           str(output_table['num_refalt'][r]),
           str(output_table['num_nonmatching'][r]),
           str(output_table['num_refalt'][r] + output_table['num_nonmatching'][r])))
    if str(output_table['category'][r]) == 'FP':
        falsepos_dedup += 1
    if str(output_table['category'][r]) == 'TP':
        truepos_dedup += 1

# print FP and TP count
# print('False positives: %s\nTrue positives: %s' % (falsepos, truepos))
print('False positives: %s\nTrue positives: %s' % (falsepos_dedup, truepos_dedup))
