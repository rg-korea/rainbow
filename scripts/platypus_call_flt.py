# Created: August 5th 2015
# Last update: July 1st 2016
# Author: Seongmin Choi <seongmin.choi@raregenomics.org>

##############
##  MODULE  ##
##############

def proc_line_flt(line, var_info, pat_id):
    # Proc
    field = line.strip().split('\t')
    chrom = field[0]
    one_pos = int(field[1])
    flt = field[6]
    chrpos = "%s:%s" % (chrom, one_pos)

    ref_base, alleles, smp_genos = var_info[chrpos] 
        
    pat_geno = smp_genos[ pat_id ]

    return (ref_base, alleles, smp_genos,
            flt, pat_geno)
# fed


############
##  INIT  ##
############

idx_s = 9 # column number from where SMP_ID starts

## CUTOFFS
#indel_length_cutoff = 7
alt_cnt_cutoff = 5 #3 #5 #2
depth_cutoff = 10 #15 #10
alt_freq_cutoff = 0.15


import sys, time
from pedigree_flt import retrieve_vcf_data
from extract_patient_vcf import proc_geno_vcf_format

prog_name = sys.argv[0].split('/')[-1]

if len(sys.argv) == 3:
    in_vcf = sys.argv[1]
    pat_id = sys.argv[2]
    print >> sys.stderr, "[%s] %s run initiated." % (time.ctime(), prog_name)
else:
    sys.exit("\nUsage: python %s <in.combined.vcf> <patient_id>\n" % prog_name)
# fi

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

    tmp_gt = pat_geno[0] # skip invariant patient genotype
    if tmp_gt == "./." or tmp_gt == "0/0":
        continue
    
    (gt, gl, gof, gq, nr, nv, abgt, inh) = proc_geno_vcf_format(pat_geno) # get allele info
     
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
        alt_idx = gt[:]
        alleles = [ alleles[x-1] for x in alt_idx ]
        depth_cnt =  [nr[x-1] for x in alt_idx]
        allele_cnt = [nv[x-1] for x in alt_idx]
        depth = max(depth_cnt)
        alt_cnt = max(allele_cnt)
    # fi
        
    min_alt_size = min([len(x) for x in alleles])

    ## FILTER
    if flt not in ["PASS"]:
        continue

#    # Indel size filter
#    indel_size = abs(len(ref_base)-min_alt_size)
#    if indel_size > indel_length_cutoff:
#        continue

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

