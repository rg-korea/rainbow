# From input VCF, select variants that corresponds to an input pedigree.
# Pedigree examples: de_novo, recessive, dominant, X-linked, Y-linked, 
# mitochondrial, compound_heterozygote. Only the former two pedigrees 
# Last update: July 1st 2016
# Author: Seongmin Choi <seongmin.choi@raregenomics.org>

###############
##  GLOBALS  ##
###############

idx_s = 9 # column number from where sample info (e.g. lca_IV-1) starts

max_alt_other_for_de_novo = 0 # max alt_cnt allowed for the patient_var_pos in parents
max_alt_freq_other_for_de_novo = 0.02 # max alt_freq allowed for the patient_var_pos in parents where the depth of parent_alt_pos is over min_depth_other_for_alt_de_novo
min_depth_other_for_alt_de_novo = 50 # cutoff for high-depth alt_pos in parents, to determine if the alt_freq in parents are low enough
min_depth_other_for_de_novo = 10 # min_depth for parent_alt_pos
de_novo_flts = (max_alt_other_for_de_novo, 
                max_alt_freq_other_for_de_novo, 
                min_depth_other_for_de_novo,
                min_depth_other_for_alt_de_novo) # tuple of filter cutoffs


###############
##  MODULES  ##
###############

def retrieve_vcf_data(vcf_file, idx_s):
    """Process VCF data and the allele genotype information
    for all input samples."""
    var_info = dict()
    for line in open(vcf_file, "r"):
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            smp_ids = line.strip().split('\t') [idx_s:] # get sample names
            continue
        else:
            field = line.strip().split('\t')
            chrom = field[0] # chromosome
            one_pos = int(field[1]) # one-based position
            chrpos = "%s:%s" % (chrom, one_pos) # chr:pos allele position
            ref_base = field[3] # C; reference base
            alt_base = field[4] # T,AC ; alterated alleles
            alleles = alt_base.split(',') # [C,T,AG]; alterated alleles
            
            formats = [ x.split(':') for x in field[idx_s:] ] # genotype information of chr:pos, for all samples
            smp_genos = { x: formats[smp_ids.index(x)] for x in smp_ids } # smp_genos[sample_id] = allele_genotype_information
            smp_genos["-"] = "0/1/2/3/4/5/6:0.0,0.0,0.0:0:0:0:0".split(':') # if sample not present
            
            var_tuple = ( ref_base, alleles, smp_genos )
            var_info[chrpos] = var_tuple # var_info[chr_pos] = variant information for the chr:pos
        # fi
    # for line end
    return var_info 
# fed

def is_pat_alt_not_in_indiv( pat_geno, indiv_geno ):
    """Determines if the patient's allele genotype information
    does not have its origin from the given individual's allele genotype information."""
    #pgt, pgl, pgof, pgq, pdepth, pnalt = pat_geno # init
    pgt, pgl, pgof, pgq, pdepth, pnalt, pagt = pat_geno # init
    #gt, gl, gof, gq, depth, nalt = indiv_geno
    gt, gl, gof, gq, depth, nalt, agt = indiv_geno

    not_in_indiv = False

    if gt == './.': # if indiv geno non-det
        return False

    pat_alleles = pgt.split('/') # patient alleles
    indiv_alleles = gt.split('/') # indiv alleles
    
    for pat_allele in pat_alleles:
        if pat_allele == '0': # ignore patient reference allele
            continue
        elif pat_allele in indiv_alleles: # if pat allele in indiv allele
            return False # even 1 pat_al present in indiv_al
        else:
            not_in_indiv = True
        # fi
    # for end

    return not_in_indiv

#    if gt == './.': # if indiv geno non-det
#        return False
#
#    pat_alleles = pgt.split('/') # patient alleles
#    indiv_alleles = gt.split('/') # indiv alleles
#    
#    for pat_allele in pat_alleles:
#        if pat_allele == '0': # ignore patient reference allele
#            continue
#        elif pat_allele in indiv_alleles: # if pat allele in indiv allele
#            return False
#        else:
#            return True
#        # fi
#    # for end
#
#    sys.exit("ERROR: patient allele: %s, indiv allele: %s"
#        % (pgt, gt))
# fed

def is_pat_alt_in_indiv( pat_geno, indiv_geno ):
    """Determines if the patient's allele genotype information
    has its origin from the given individual's allele genotype information."""
    pgt, pgl, pgof, pgq, pdepth, pnalt = pat_geno # init
    gt, gl, gof, gq, depth, nalt = indiv_geno

    if gt == './.': # if indiv geno non-det
        return False

    pat_alleles = pgt.split('/') # patient alleles
    indiv_alleles = gt.split('/') # indiv alleles
    
    for pat_allele in pat_alleles:
        if pat_allele == '0': # ignore patient reference allele
            continue
        elif pat_allele in indiv_alleles: # if pat allele in indiv allele
            return True
        else:
            return False
        # fi
    # for end

    sys.exit("ERROR: patient allele: %s, indiv allele: %s"
        % (pgt, gt))
# fed

def is_indiv_alt_cnt_infrequent( indiv_geno, de_novo_flts ):
    """Determines if the individual's allele genotype information
    shows a low frequency of variants."""
    #gt, gl, gof, gq, nr, nv = indiv_geno
    gt, gl, gof, gq, nr, nv, agt = indiv_geno

    (max_alt_other_for_de_novo, 
     max_alt_freq_other_for_de_novo, 
     min_depth_other_for_de_novo,
     min_depth_other_for_alt_de_novo) = de_novo_flts

    if gt == './.': # if indiv geno non-det
        return False

    gt = [int(x) for x in gt.split('/')] # indiv alleles
    nr = [int(x) for x in nr.split(',')]
    nv = [int(x) for x in nv.split(',')]

    max_depth = max(nr)
    if max_depth < min_depth_other_for_de_novo:
        return False

    min_alt_cnt = min(nv)
    alt_freq = float(min_alt_cnt) / max_depth

    if min_alt_cnt > max_alt_other_for_de_novo:
        if (alt_freq < max_alt_freq_other_for_de_novo and
            max_depth >= min_depth_other_for_alt_de_novo):
            return True
        else:
            return False
        # fi
    else:
        return True
# fed

def is_indiv_geno_1_alt_hetero( indiv_geno ):
    """Determines if the individual's allele genotype information
    has two heterozygotic alterations; e.g. ref:A/A, indiv:C/T"""
    gt, gl, gof, gq, depth, nalt = indiv_geno

    if gt == './.': # if indiv geno non-det
        return False

    if ('0' in gt) and gt.count('0') == 1:
        return True
    else:
        return False

    sys.exit("ERROR: indiv allele: %s" % gt)
# fed

def is_indiv_geno_2_alt_hetero( indiv_geno ):
    """Determines if the individual's allele genotype information
    has two heterozygotic alterations; e.g. ref:A/A, indiv:C/T"""
    gt, gl, gof, gq, depth, nalt = indiv_geno

    if gt == './.': # if indiv geno non-det
        return False

    if '0' in gt:
        return False

    set_alleles = set( gt.split('/') ) # set of indiv alleles
    if len(set_alleles) == 1:
        return False
    else:
        return True
    # fi

    sys.exit("ERROR: indiv allele: %s" % gt)
# fed

def is_indiv_geno_alt_homo( indiv_geno ):
    """Determines if the individual's allele genotype information
    has a homozygotic alteration; e.g. ref:A/A, indiv:T/T"""
    gt, gl, gof, gq, depth, nalt = indiv_geno

    if gt == './.': # if indiv geno non-det
        return False

    if '0' in gt:
        return False

    set_alleles = set( gt.split('/') ) # set of indiv alleles
    if len(set_alleles) == 1:
        return True
    else:
        return False
    # fi

    sys.exit("ERROR: indiv allele: %s" % gt)
# fed

def is_indiv_geno_var( indiv_geno ):
    """Determines if the individual's allele genotype information
    retains a variant."""
    gt, gl, gof, gq, depth, nalt = indiv_geno
    if gt in ['./.', '0/0']:
        return False
    elif gt.count('.'):
        sys.exit("ERROR: GT = %s" % gt)
    else:
        return True
    # fi
# fed

def is_alts_from_mot_and_fat(pat_gt, fat_gt, mot_gt):
    """Determines if the patient's allele genotype information
    has its origins from both mother and father."""
    pat_alleles = pat_gt.split('/') # init
    fat_alleles = fat_gt.split('/')
    mot_alleles = mot_gt.split('/')

    is_in_fat, is_in_mot = False, False # init

    for i, pat_allele in enumerate(pat_alleles): # is in mot and fat?
        if pat_allele in fat_alleles: # found in fat?
            is_in_fat = True
            fat_idx = fat_alleles.index(pat_allele)
            pat_alleles[i] = None
            fat_alleles[fat_idx] = None
            continue
        # fi
        if pat_allele in mot_alleles: # found in mot?
            is_in_mot = True
            mot_idx = mot_alleles.index(pat_allele)
            pat_alleles[i] = None
            mot_alleles[mot_idx] = None
            continue
        # fi
    # for end

    return is_in_fat, is_in_mot
# fed

def is_indiv_alt_recessive( pat_geno, fat_geno, mot_geno ):
    # Proc
    pat_gt, pat_gl, pat_gof, pat_gq, pat_depth, pat_nalt = pat_geno
    fat_gt, fat_gl, fat_gof, fat_gq, fat_depth, fat_nalt = fat_geno
    mot_gt, mot_gl, mot_gof, mot_gq, mot_depth, mot_nalt = mot_geno

    # Only variants in patient
    is_pat_geno_var = is_indiv_geno_var( pat_geno )
    if not is_pat_geno_var:
        return False

    # Is pat alt hetero or homo
    is_pat_var_hetero = is_indiv_geno_2_alt_hetero( pat_geno ) # 1/2?
    is_pat_var_homo = is_indiv_geno_alt_homo( pat_geno ) # 1/1?
    if not (is_pat_var_hetero or is_pat_var_homo):
        return False

    # If pat genotype same as either mot or fat
    if pat_gt in [fat_gt, mot_gt]:
        return False

    # Test segregation: recessive
    is_in_fat, is_in_mot = is_alts_from_mot_and_fat(pat_gt, fat_gt, mot_gt)
    is_segregated = (is_in_fat and is_in_mot)

    return is_segregated
# fed

def print_de_novo(in_vcf, indiv_ids, idx_s, de_novo_flts):
    # Make prior
    var_info = retrieve_vcf_data(in_vcf, idx_s)

    # Proc VCF
    for line in open(in_vcf):
        # Print '#' lines
        if line.startswith("#"):
            print line.strip()
            continue

        # Proc line
        (ref_base, alleles, smp_genos,
         pat_geno, fat_geno, mot_geno) = proc_line(line, var_info, indiv_ids)

        # Only variants in patient
        if not is_indiv_geno_var( pat_geno ):
            continue
        
        # If patient variant not in father or mother, 
        not_in_fat = is_pat_alt_not_in_indiv( pat_geno, fat_geno )
        not_in_mot = is_pat_alt_not_in_indiv( pat_geno, mot_geno )

        # Select only parents with infrequent variant alleles
        no_alt_in_fat = is_indiv_alt_cnt_infrequent( fat_geno, de_novo_flts )
        no_alt_in_mot = is_indiv_alt_cnt_infrequent( mot_geno, de_novo_flts )

        is_de_novo = (( not_in_fat and not_in_mot ) and # is var de novo?
                      ( no_alt_in_fat and no_alt_in_mot ))

        #is_de_novo = ( not_in_fat and not_in_mot )# is var de novo?

        # Print line if de novo
        if is_de_novo:
            print line.strip()
    # for end
# fed

def print_recessive(in_vcf, indiv_ids, idx_s):
    # Make prior
    var_info = retrieve_vcf_data(in_vcf, idx_s)

    # Open VCF
    for line in open(in_vcf):
        # Print '#' lines
        if line.startswith("#"):
            print line.strip()
            continue

        # Proc line
        (ref_base, alleles, smp_genos,
         pat_geno, fat_geno, mot_geno) = proc_line(line, var_info, indiv_ids)

        # If patient var is in father and mother
        is_recessive = is_indiv_alt_recessive( pat_geno, fat_geno, mot_geno )

        # Print line if recessive
        if is_recessive:
            print line.strip()
    # for end
# fed

 
############
##  MAIN  ##
############

def get_abstract_gt( indiv_geno ):
    gt, gl, gof, gq, nr, nv = indiv_geno
    if gt == "./.": # not defined
        agt = "-"
    elif gt == "0/0": # reference
        agt = "ref"
    elif is_indiv_geno_1_alt_hetero( indiv_geno ):
        agt = "het1"
    elif is_indiv_geno_2_alt_hetero( indiv_geno ):
        agt = "het2"
    elif is_indiv_geno_alt_homo( indiv_geno ):
        agt = "hom"
    elif len(gt) == 13:
        agt = "parent" # gene
    else:
        sys.exit("ERROR: indiv_geno = %s" % indiv_geno)
    # fi

    return agt
# fed

def get_inh_status( pat_geno, fat_geno, mot_geno, de_novo_flts ):
    pat_gt, pat_gl, pat_gof, pat_gq, pat_nr, pat_nv, pat_agt = pat_geno # unpack
    fat_gt, fat_gl, fat_gof, fat_gq, fat_nr, fat_nv, fat_agt = fat_geno
    mot_gt, mot_gl, mot_gof, mot_gq, mot_nr, mot_nv, mot_agt = mot_geno

    field = line.strip().split('\t')
    
    if fat_gt == "./." or mot_gt == "./.": # either parent gt is not determined
        return "-" # return NA's
    
    if len(fat_gt) == 13 or len(mot_gt) == 13: # father sample x, mother sample o
        return "inherited" # inherited
    else:
        if pat_agt == "ref": # patient abstract gt is ref
            return "inherited" # inherited
        else: # patient abstract gt is non-ref
            not_in_fat = is_pat_alt_not_in_indiv( pat_geno, fat_geno )
            not_in_mot = is_pat_alt_not_in_indiv( pat_geno, mot_geno )
            no_alt_in_fat = is_indiv_alt_cnt_infrequent( fat_geno, de_novo_flts )
            no_alt_in_mot = is_indiv_alt_cnt_infrequent( mot_geno, de_novo_flts )
            is_de_novo = (( not_in_fat and not_in_mot ) and # is var de novo?
                          ( no_alt_in_fat and no_alt_in_mot ))

            if is_de_novo:
                return "de_novo" # de novo
            else:
                return "inherited" # inherited
            # fi
        # fi
    # fi

    sys.exit("ERROR: pat_geno:%s, fat_geno:%s, mot_geno:%s" 
        % (pat_geno, fat_geno, mot_geno))
# fed

def proc_line(line, var_info, indiv_ids, de_novo_flts):
    # Proc
    field = line.strip().split('\t')
    chrom = field[0]
    one_pos = int(field[1])
    chrpos = "%s:%s" % (chrom, one_pos)

    ref_base, alleles, smp_genos = var_info[chrpos] 
       
    pat_ids, fat_id, mot_id = indiv_ids # sample_id or na
    pat_genos = [smp_genos[x] for x in pat_ids] # get patient allele info; e.g. 0/1:-20.37,0.0,-37.37:35:99:22:5
    fat_geno = smp_genos[ fat_id ] # get father allele info
    mot_geno = smp_genos[ mot_id ] # get mother allele info

    pat_agts = [get_abstract_gt(x) for x in pat_genos]
    fat_agt = get_abstract_gt( fat_geno )
    mot_agt = get_abstract_gt( mot_geno )

    for i, pat_agt in enumerate(pat_agts):
        pat_genos[i] += [ pat_agt ]
    fat_geno += [ fat_agt ]
    if mot_agt not in mot_geno:
        mot_geno += [ mot_agt ]

    pat_inhs = [get_inh_status(x, fat_geno, mot_geno, de_novo_flts) for x in pat_genos]

    for i, pat_inh in enumerate(pat_inhs):
        pat_genos[i] += [ pat_inh ]

    return (ref_base, alleles,
            pat_genos, fat_geno, mot_geno)
# fed


if __name__ == "__main__":

    # Get system arguments
    import sys, time
    prog_name = sys.argv[0].split('/')[-1]
    if len(sys.argv) == 4:
        in_vcf = sys.argv[1]
        fat_id = sys.argv[2]
        mot_id = sys.argv[3]
        print >> sys.stderr, "[%s] %s run initiated." % (time.ctime(), prog_name)
    else:
        sys.exit(("\nUsage: python %s <in.vcf> <fat_id> <mot_id>\n" % prog_name))
    # fi

    #print_vcf( in_vcf, pat_info, fat_info, mot_info, idx_s )i

    var_info = retrieve_vcf_data(in_vcf, idx_s)
    for line in open(in_vcf):
        if line.startswith("##"):
            print line.strip()
            continue
        elif line.startswith('#'):
            field = line.strip().split('\t')
            pat_ids = [x for x in field[idx_s:] if (x!=fat_id and x!=mot_id)]
            fat_idx = field.index( fat_id ) if fat_id != "-" else None
            mot_idx = field.index( mot_id ) if mot_id != "-" else None
            pat_idxs = [field.index(x) for x in pat_ids]
            indiv_ids = (pat_ids, fat_id, mot_id)
            print "##FORMAT=<ID=AG,Number=.,Type=String,Description=\"Abstract genotype: -, ref, het1, het2, hom, not_ref\">"
            print "##FORMAT=<ID=IH,Number=.,Type=String,Description=\"Inheritence status: inherited, de_novo\">"
            print line.strip()
            continue

        # Proc line
        (ref_base, alleles,
         pat_geno, fat_geno, mot_geno) = proc_line(line, var_info, indiv_ids, de_novo_flts)

        field = line.strip().split('\t')
        field[idx_s-1] = "GT:GL:GOF:GQ:NR:NV:AG:IH"
        for i, pat_idx in enumerate(pat_idxs):
            field[pat_idx] = ':'.join(pat_geno[i])
        if fat_idx:
            field[fat_idx] = ':'.join(fat_geno)
        if mot_idx:
            field[mot_idx] = ':'.join(mot_geno)
        print '\t'.join(field)
    # for line end



