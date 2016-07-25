# Created: August 5th 2015
# Last update: July 1st 2016
# Author: Seongmin Choi <seongmin.choi@raregenomics.org>

############
##  INIT  ##
############

idx_s = 9 # column number from where SMP_ID starts

## CUTOFFS
alt_cnt_cutoff = 5 #3 #5 #2
depth_cutoff = 10 #15 #10
alt_freq_cutoff = 0.15


import sys, time
from specify_variants import retrieve_vcf_data
from extract_patient_vcf import proc_geno_vcf_format

prog_name = sys.argv[0].split('/')[-1]

if len(sys.argv) == 3:
    in_vcf = sys.argv[1]
    pat_id = sys.argv[2]
    print >> sys.stderr, "[%s] %s run initiated." % (time.ctime(), prog_name)
else:
    sys.exit("\nUsage: python %s <in.combined.vcf> <patient_id>\n" % prog_name)
# fi

##############
##  MODULE  ##
##############

def proc_line_flt(line, var_info, pat_id):
    # Proc
    field = line.strip().split('\t')
    chrom = field[0]
    one_pos = int(field[1])
    ref_base = field[3]
    alt_base = field[4]
    flt = field[6]
    var_key = "%s:%s:%s:%s" % (chrom, one_pos, ref_base, alt_base)

    ref_base, alleles, smp_genos = var_info[var_key] 
        
    pat_geno = smp_genos[ pat_id ]

    return (ref_base, alleles, smp_genos,
            flt, pat_geno)
# fed

##@## def get_depth_and_alt_cnt():

############
##  MAIN  ##
############

# Make prior
var_info = retrieve_vcf_data(in_vcf, idx_s)

for line in open(in_vcf, "r"):
    if line.startswith('#'):
        print line.strip()
        continue

    # Proc line
    (ref_base, alleles, smp_genos, flt,
    pat_geno) = proc_line_flt(line, var_info, pat_id)
    if pat_geno["GT"] in ["./.", "0/0"]: # filter out non-variants
        continue

    (gt, nr, nv, abgt, inh) = proc_geno_vcf_format(pat_geno) # get allele info
    # gt: [0,1] / [1,2] / [0,2] / ...
    # nr: [45,45] (depth info per allele)
    # nv: [16,29]: 16 reads supporting gt[0] and 29 reads supporting gt[1]
     
    # Get depth and altered allele counts
    if 0 in gt: # patient is heterozygous 0/1
        alt_idx = filter(lambda x: x != 0, gt) [0]
        alleles = [ alleles[alt_idx-1] ]
        depth = nr[alt_idx-1]
        alt_cnt = nv[alt_idx-1]
    elif len(set(gt)) == 1: # patient is homozygous
        alt_idx = gt[0]
        alleles = [ alleles[alt_idx-1] ]
        depth = nr[alt_idx-1]
        alt_cnt = nv[alt_idx-1]
    else: # patient is heterozygous 1/2
        alt_idx = gt[:] # two alt idxes: [1,2] (or [2,3], [3,4] ...)
        alleles = [ alleles[x-1] for x in alt_idx ]
        depth_cnt =  [nr[x-1] for x in alt_idx]
        allele_cnt = [nv[x-1] for x in alt_idx]
        depth = max(depth_cnt)
        alt_cnt = max(allele_cnt)
    # fi
        
    min_alt_size = min([len(x) for x in alleles])

    ## FILTERS
    # "FILTER" column filter
    if flt not in ["PASS"]:
        continue

    # Variant count filter
    if alt_cnt < alt_cnt_cutoff:
        continue

    # Depth cutoff
    if depth < depth_cutoff:
        continue

    # Variant frequency cutoff:
    alt_freq = alt_cnt / float(depth)
    if alt_freq < alt_freq_cutoff:
        continue

    # Print line
    print line.strip()
    # fi
# for line end

