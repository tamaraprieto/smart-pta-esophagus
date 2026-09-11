#!/usr/bin/env python
import argparse
import sys
import pysam
import pandas
import os   

parser = argparse.ArgumentParser(description=' ')
parser.add_argument('--bamfile')
parser.add_argument('--bedfile')
parser.add_argument('--min_mapping_q', type=int, default=20)
parser.add_argument('--min_base_q', type=int, default=15)
parser.add_argument('--max_edit_dist', type=int, default=5)
args = parser.parse_args()

bamfile = args.bamfile
bedfile = args.bedfile
ref_fasta = pysam.FastaFile('/gpfs/commons/groups/landau_lab/tprieto/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta')

min_mapping_q = args.min_mapping_q # default 20
min_base_q = args.min_base_q # default 15 
max_edit_dist = args.max_edit_dist # default 5 


# check that expanded_hetsnps bedfile is not empty, if empty, exit
if os.stat(bedfile).st_size == 0:
    print('Bedfile is empty')
    sys.exit()
    
# generate table to store the possible base combos present in the reads of interest
# first column is the base present in position 1, second column is the base present in position 2
# third column is the number of reads containing that combo
rows, cols = (16, 3)
basecount = [[0 for i in range(cols)] for j in range(rows)]
bases = ['A', 'T', 'C', 'G']
index = 0
for i in bases:
    for j in bases:
        basecount[index][0] = i
        basecount[index][1] = j
        index+=1
    
    
# create table from bed file that stores position of each variant and the position of the corresponding hetSNP
bed=pandas.read_table(bedfile,sep='\t',header=None)
positions = pandas.DataFrame(list(zip(bed.iloc[:,0],bed.iloc[:,1],bed.iloc[:,5],bed.iloc[:,6],bed.iloc[:,8],bed.iloc[:,10],bed.iloc[:,11])), columns=['chr', 'variant', 'snv_ref', 'snv_alt', 'hetSNP', 'hetSNP_ref', 'hetSNP_alt'])
positions['variant'] += 100 # we expanded the bed file intervals by 100 so +100 to get original pos

# arrays for counting reads that support each combo of variant/hetSNP
hetsnp_ref_reads = [[0 for i in range(3)] for j in range(2)]
hetsnp_alt_reads = [[0 for i in range(3)] for j in range(2)]


# for each row in the bed, count the bases at the variant and the bases at the hetsnp
# then intersect based on read name
bam = pysam.AlignmentFile(bamfile, "rb")

# iterate over our list of snvs
for r in positions.index:
    chr = positions['chr'][r]
    varposition = positions['variant'][r]
    hetsnpposition = positions['hetSNP'][r]
    
    # store the variant ref/alt and hetSNP ref/alt
    hetsnp_ref = positions['hetSNP_ref'][r]
    hetsnp_alt = positions['hetSNP_alt'][r]
    snv_ref = positions['snv_ref'][r]
    snv_alt = positions['snv_alt'][r]
    # fill in the arrays with the bases
    hetsnp_ref_reads[0][0] = hetsnp_ref
    hetsnp_ref_reads[1][0] = hetsnp_ref
    hetsnp_ref_reads[0][1] = snv_ref
    hetsnp_ref_reads[1][1] = snv_alt
    hetsnp_ref_reads[0][2] = 0
    hetsnp_ref_reads[1][2] = 0
    hetsnp_alt_reads[0][0] = hetsnp_alt
    hetsnp_alt_reads[1][0] = hetsnp_alt
    hetsnp_alt_reads[0][1] = snv_ref
    hetsnp_alt_reads[1][1] = snv_alt
    hetsnp_alt_reads[0][2] = 0
    hetsnp_alt_reads[1][2] = 0
    other_reads = basecount  # for reads not matching the ref/alt
    
    # create table to store only reads which cover both bases
    linkedreads = pandas.DataFrame(columns=['readname', 'variant_base', 'hetsnp_base'])
    
    # get pileup of reads covering variant
    for pileupcolumn in bam.pileup(chr, varposition, min_mapping_quality = min_mapping_q, min_base_quality = min_base_q): # filter out low mapping quality and low base quality (at variant) (default = 13)
        if pileupcolumn.pos == varposition:
            print("\ncoverage at position %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.get_tag("NM") < max_edit_dist: # filter out edit distance > 5 (NM tag) 
                    insert_row = {
                        "readname": pileupread.alignment.query_name,
                        "variant_base": pileupread.alignment.query_sequence[pileupread.query_position],
                        "hetsnp_base": None
                    }
                    linkedreads = pandas.concat([linkedreads, pandas.DataFrame([insert_row])])
                    # print(pileupread.query_position, pileupread.alignment.query_sequence[pileupread.query_position])
    
    # debugging: 
    # print("het snp: ",  hetsnpposition, "\tmapping q: ", min_mapping_q, "\tmin base q: ", min_base_q)
    
    # get pileup of reads covering hetsnp
    for pileupcolumn in bam.pileup(chr, hetsnpposition, min_mapping_quality = min_mapping_q): # filter out low mapping quality and low base quality (at variant) (default = 13)
        if pileupcolumn.pos == hetsnpposition - 1:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.get_tag("NM") < max_edit_dist: # filter out edit distance > 5 (NM tag) 
                    linkedreads.loc[linkedreads['readname'] == pileupread.alignment.query_name, 'hetsnp_base'] = pileupread.alignment.query_sequence[pileupread.query_position]
                    # print(pileupread.query_position, pileupread.alignment.query_sequence[pileupread.query_position])
                    
                    # if the read covers both positions, print the info from the row 
                    # if linkedreads['readname'] == pileupread.alignment.query_name:
                    if linkedreads['readname'].str.contains(pileupread.alignment.query_name).any():
                        variant_base = linkedreads.loc[linkedreads['readname'] == pileupread.alignment.query_name, 'variant_base'].values[0]
                        hetsnp_base = linkedreads.loc[linkedreads['readname'] == pileupread.alignment.query_name, 'hetsnp_base'].values[0]
                        print('\tbase in %s, read %s at position %s = %s, position %s = %s' % 
                            (chr, pileupread.alignment.query_name, 
                            varposition, 
                            # linkedreads.loc[linkedreads['readname'] == pileupread.alignment.query_name, 'variant_base'].values[0],
                            variant_base,
                            hetsnpposition, 
                            # linkedreads.loc[linkedreads['readname'] == pileupread.alignment.query_name, 'hetsnp_base'].values[0]))
                            hetsnp_base))
                        
                        # if hetsnp has ref allele
                        if hetsnp_base == hetsnp_ref:
                            if variant_base == snv_ref:
                                hetsnp_ref_reads[0][2] += 1
                            elif variant_base == snv_alt:
                                hetsnp_ref_reads[1][2] += 1
                            else:
                                # add to basecount
                                for i in range(rows):
                                    if str(basecount[i][0]) + str(basecount[i][1]) == variant_base + hetsnp_base:
                                        basecount[i][2] += 1
                        # if hetsnp has alt allele
                        elif hetsnp_base == hetsnp_alt:
                            if variant_base == snv_ref:
                                hetsnp_alt_reads[0][2] += 1
                            elif variant_base == snv_alt:
                                hetsnp_alt_reads[1][2] += 1
                            else:
                                # add to basecount
                                for i in range(rows):
                                    if str(basecount[i][0]) + str(basecount[i][1]) == variant_base + hetsnp_base:
                                        basecount[i][2] += 1
                        
                        else:
                            # add to basecount 
                            for i in range(rows):
                                if str(basecount[i][0]) + str(basecount[i][1]) == variant_base + hetsnp_base:
                                    basecount[i][2] += 1
    # print(linkedreads)
    linkedreads = linkedreads.dropna()
    
    # print total counts for current variant
    # while also clearing counts before iterating onto next variant
    print('Variant at position %s' % varposition)
    
    if hetsnp_ref_reads[0][2] > 0:
        print('\tReads supporting REF %s at variant and REF %s at hetSNP: %s' % (str(hetsnp_ref_reads[0][1]), str(hetsnp_ref_reads[0][0]), str(hetsnp_ref_reads[0][2])))
    hetsnp_ref_reads[0][2] = 0
    if hetsnp_ref_reads[1][2] > 0:
        print('\tReads supporting ALT %s at variant and REF %s at hetSNP: %s' % (str(hetsnp_ref_reads[1][1]), str(hetsnp_ref_reads[1][0]), str(hetsnp_ref_reads[1][2])))
    hetsnp_ref_reads[1][2] = 0
    
    # print counts for alt hetsnp
    if hetsnp_alt_reads[0][2] > 0:
        print('\tReads supporting REF %s at variant and ALT %s at hetSNP: %s' % (str(hetsnp_alt_reads[0][1]), str(hetsnp_alt_reads[0][0]), str(hetsnp_alt_reads[0][2])))
    hetsnp_alt_reads[0][2] = 0
    if hetsnp_alt_reads[1][2] > 0:
        print('\tReads supporting ALT %s at variant and ALT %s at hetSNP: %s' % (str(hetsnp_alt_reads[1][1]), str(hetsnp_alt_reads[1][0]), str(hetsnp_alt_reads[1][2])))
    hetsnp_alt_reads[1][2] = 0
    
    # print counts for reads not matching ref/alt
    print('\tNon-matching reads:')
    for r in range(rows):
        if basecount[r][2] > 0:
            print('\tReads with %s at variant and %s at local hetSNP: %s' % (str(basecount[r][0]), str(basecount[r][1]), str(basecount[r][2])))
        basecount[r][2] = 0

bam.close()
# print(positions)
# print(linkedreads)
