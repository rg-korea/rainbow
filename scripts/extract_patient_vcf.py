# Last update: July 1st 2016
# Author: Seongmin Choi <seongmin.choi@raregenomics.org>

############
##  INIT  ##
############

import sys

idx_s = 9 # column number from where SMP_ID starts


###############
##  MODULES  ##
###############

def make_print_line(line, print_alt, print_geno, idx_s):
    field = line.strip().split('\t')
    print_field = field[:idx_s]
    print_field[4] = print_alt
    print_field.append( print_geno )
    print_line = "\t".join(print_field)

    return print_line
# fed

def proc_geno_vcf_format(geno):
    gt = [int(x) for x in geno[0].split('/')] # [2,3]
    gl = geno[1]
    gof = geno[2]
    gq = geno[3]
    nr = [int(x) for x in geno[4].split(',')] # depths
    nv = [int(x) for x in geno[5].split(',')] # navrs
    abgt = geno[6]
    inh = geno[7]

    return (gt, gl, gof, gq, nr, nv, abgt, inh)
# fed

def get_print_nr_nv(gt, nr, nv):
    if gt == "./.":
        sys.exit("ERROR: gt = %s" % gt)
    elif gt.count('0'):
        flag_one_alt = True
        print_nr = str(nr[0])
        print_nv = str(nv[0])
    elif len(set(gt)) == 1:
        flag_one_alt = True
        print_nr = str(nr[ gt[0]-1 ])
        print_nv = str(nv[ gt[0]-1 ])
    else:
        flag_one_alt = False
        print_nr = ",".join([str(nr[x-1]) for x in gt])
        print_nv = ",".join([str(nv[x-1]) for x in gt])
    # fi
    
    return print_nr, print_nv
# fed

def get_print_geno_alleles(ref_base, alleles, geno):
    # Proc
    ref_alt = [ref_base] + alleles # [0A, 1C,2AC,3G,4T]
    (gt, gl, gof, gq, nr, nv, abgt, inh) = proc_geno_vcf_format(geno)

    print_nr, print_nv = get_print_nr_nv(gt, nr, nv) # NR and NV to print

    gt_based_alleles = [ref_alt[x] for x in gt]
    print_gt = gt[:] 
    print_alt = [ a for a in alleles if a in gt_based_alleles ] # ALT to print

    # GT to print
    for input_idx, gt_based_allele in enumerate(gt_based_alleles):
        start_idx = ref_alt.index( gt_based_allele )
        print_gt_idx = start_idx
        for i in xrange(start_idx-1, 0, -1):
            if ref_alt[i] in gt_based_alleles:
                continue
            else:
                print_gt_idx -= 1
            # fi
        # for end
        print_gt[input_idx] = print_gt_idx
    # for end

    print_gt = "/".join( [str(x) for x in print_gt] )
    print_geno = ':'.join([print_gt, gl, gof, gq, print_nr, print_nv, abgt, inh])
    print_alt = ",".join( print_alt )

    return print_geno, print_alt
# fed

def proc_line_regard_patient(line, var_info, pat_id):
    # Proc
    field = line.strip().split('\t')
    chrom = field[0]
    one_pos = int(field[1])
    chrpos = "%s:%s" % (chrom, one_pos)

    ref_base, alleles, smp_genos = var_info[chrpos] 
        
    pat_geno = smp_genos[ pat_id ]

    return (ref_base, alleles, smp_genos,
            pat_geno)
# fed


############
##  MAIN  ##
############

if __name__ == "__main__":

    # Get system arguments
    import sys, time
    from pedigree_flt import retrieve_vcf_data
    prog_name = sys.argv[0].split('/')[-1]
    if len(sys.argv) == 3:
        in_vcf = sys.argv[1]
        pat_id = sys.argv[2] # patient ID
        print >> sys.stderr, "[%s] %s run initiated." % (time.ctime(), prog_name)
    else:
        sys.exit(("\nUsage: python %s <in.combined.vcf> <patient_id>\n" % prog_name))
    # fi
    
    # Make priors
    var_info = retrieve_vcf_data(in_vcf, idx_s)
    
    # Proc VCF
    for line in open(in_vcf):
        # Print '#' lines
        if line.startswith("##"):
            print line.strip()
            continue
        elif line.startswith("#CHROM"):
            field = line.strip().split('\t')
            to_print = field[:idx_s] + [pat_id]
            print "\t".join(to_print)
            continue
        
        # Proc line
        (ref_base, alleles, smp_genos, 
        pat_geno) = proc_line_regard_patient(line, var_info, pat_id)
        
        # Print patient VCF line output
        print_geno, print_alt = get_print_geno_alleles(ref_base, alleles, pat_geno)
        print_line = make_print_line(line, print_alt, print_geno, idx_s)
        print print_line
    # for end
 
