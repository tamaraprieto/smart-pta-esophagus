import pysam
import argparse
import os

# ----------------------
# Parse arguments
# ----------------------
parser = argparse.ArgumentParser(description="Count total and alternative allele reads per SNV in a BAM file.")
parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
parser.add_argument("-r", "--ref", required=True, help="Reference FASTA file")
parser.add_argument("-s", "--snv", required=True, help="SNV list file (CHROM START POS REF ALT)")
parser.add_argument("-o", "--outdir", default=".", help="Output directory (default: current)")
args = parser.parse_args()

# ----------------------
# Set variables
# ----------------------
bam_file = args.bam
ref_fasta = args.ref
snv_file = args.snv
outdir = args.outdir

# Derive sample name from BAM
sample_name = os.path.basename(bam_file).replace(".bam", "")

# Output files
depth_file = os.path.join(outdir, f"{sample_name}_total_depth.tsv")
alt_file = os.path.join(outdir, f"{sample_name}_alt_count.tsv")

# ----------------------
# Open BAM and reference
# ----------------------
bam = pysam.AlignmentFile(bam_file, "rb")
ref = pysam.FastaFile(ref_fasta)

# ----------------------
# Write output files
# ----------------------


with open(depth_file, "w") as depth_out, open(alt_file, "w") as alt_out:
    depth_out.write(f"Variant\t{sample_name}\n")
    alt_out.write(f"Variant\t{sample_name}\n")

    with open(snv_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            # BED format: chrom start end ref alt
            chrom, start, end, ref_base, alt_base = line.strip().split()[:5]
            pos = int(end) 

            variant_id = f"{chrom}:{pos}_{ref_base}>{alt_base}"
            #print(variant_id)		

            total_depth = 0
            alt_count = 0

            # parameters you can tune
            MIN_MAPQ = 20
            MIN_BASEQ = 15
            MAX_DEPTH = 2000000   
            FILTER_DUPLICATES = True
            FILTER_SECONDARY_SUPP = True
            #COUNT_FRAGMENTS = False  
            
            for pileupcolumn in bam.pileup(
                    chrom,
                    pos-1,
                    pos,
                    truncate=True,
                    max_depth=MAX_DEPTH,
                    min_base_quality=0):  
                if pileupcolumn.reference_pos + 1 != pos:
                    continue
            
                total_depth = 0
                alt_count = 0
            
                #seen_fragments = set()   # optional, for fragment-aware counting
            
                for pr in pileupcolumn.pileups:
                    aln = pr.alignment
            
                    # skip deletions / reference skips
                    if pr.is_del or pr.is_refskip:
                        continue
            
                    # skip if query_position is None (safety)
                    if pr.query_position is None:
                        continue
            
                    # basic alignment-level filters
                    if FILTER_DUPLICATES and aln.is_duplicate:
                        continue
                    if FILTER_SECONDARY_SUPP and (aln.is_secondary or aln.is_supplementary):
                        continue
                    if aln.mapping_quality < MIN_MAPQ:
                        continue
            
                    # base and base-quality
                    qpos = pr.query_position
                    base = aln.query_sequence[qpos].upper()
                    baseq = aln.query_qualities[qpos] if aln.query_qualities is not None else None
                    if (MIN_BASEQ is not None) and (baseq is not None) and (baseq < MIN_BASEQ):
                        continue
            
                    # passed all filters
                    total_depth += 1
                    if base == alt_base:
                        alt_count += 1
            
            # now write out total_depth and alt_count
            depth_out.write(f"{variant_id}\t{total_depth}\n")
            alt_out.write(f"{variant_id}\t{alt_count}\n")


bam.close()
ref.close()
print(f"Files written: {depth_file}, {alt_file}")

