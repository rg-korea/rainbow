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
            smp_genos = { "-": "0/1/2/3/4/5/6:0.0,0.0,0.0:0:0:0:0" } # if sample not present
            smp_genos = { x: formats[smp_ids.index(x)] for x in smp_ids } # smp_genos[sample_id] = allele_genotype_information
            
            var_tuple = ( ref_base, alleles, smp_genos )
            var_info[chrpos] = var_tuple # var_info[chr_pos] = variant information for the chr:pos
        # fi
    # for line end
    return var_info 
# fed

def get_abgt_list(indiv_flt):
    if indiv_flt == "not_ref":
        return ["het1", "het2", "hom"]
    elif indiv_flt == "all":
        return ["ref", "het1", "het2", "hom"]
    elif indiv_flt == "het":
        return ["het1", "het2"]
    elif indiv_flt == "not_hom":
        return ["ref", "het1", "het2"]
    elif indiv_flt == "het1_hom":
        return ["het1", "hom"]
    elif indiv_flt == "ref_hom":
        return ["het1", "hom"]
    else:
        return [indiv_flt]
    # fi
# fed

def chrom_fit_pedigree(inh_type, chrom):
    chrom_id = chrom.replace("chr",'')
    if inh_type.startswith('X'):
        if chrom_id.count('X'):
            return True
        else:
            return False
    elif inh_type.startswith('Y'):
        if chrom_id.count('Y'):
            return True
        else:
            return False
    elif inh_type.startswith('A'):
        if (chrom_id.count('X') or chrom_id.count('Y') 
            or chrom_id.count('M') or chrom_id.count("MT")):
            return False # if is X, Y, or M/MT(mitochondrial)
        else:
            return True # autosomes
    elif inh_type == "MT": # mitochondrial:
        if chrom_id.count('M') or chrom_id.count("MT"):
            return True
        else:
            return False
    else:
        sys.exit("ERROR: inh_type: %s, chromosome:%s" % (inh_type, chrom))
    # fi
# fed

def print_filtered_vcf(in_vcf, inh_type, pat_id, fat_id, mot_id,
                       pat_abgt_list, fat_abgt_list, mot_abgt_list):
    for line in open(in_vcf, "r"):
        if line.startswith("##"): # header line
            print line.strip()
            continue
        elif line.startswith('#'): # column name line
            field = line.strip().split('\t')
            pat_idx = field.index(pat_id)
            fat_idx = field.index(fat_id) if fat_id != "-" else None
            mot_idx = field.index(mot_id) if mot_id != "-" else None
            print line.strip()
            continue
        else: # variant line
            field = line.strip().split('\t')
            chrom = field[0]
            if not chrom_fit_pedigree(inh_type, chrom): # flt out non-fit chromosomes
                continue
            pat_abgt, pat_inh = field[pat_idx].split(':')[-2:]
            fat_abgt = field[fat_idx].split(':')[-1] if fat_idx else "-"
            mot_abgt = field[mot_idx].split(':')[-1] if mot_idx else "-"
            if pat_inh != inh_flt: # filter according to patient inheritance
                continue
            if pat_abgt not in pat_abgt_list:
                continue
            if fat_abgt not in fat_abgt_list:
                continue
            if mot_abgt not in mot_abgt_list:
                continue
            print line.strip()                
        #  fi
    # for line end
# fed

############
##  MAIN  ##
############

if __name__ == "__main__":

    # Get system arguments
    import sys, time
    prog_name = sys.argv[0].split('/')[-1]
    if len(sys.argv) == 7:
        in_vcf = sys.argv[1]
        inh_flt = sys.argv[2]
        inh_type = sys.argv[3]
        pat_id, pat_flt = sys.argv[4].split(':')
        fat_id, fat_flt = sys.argv[5].split(':')
        mot_id, mot_flt = sys.argv[6].split(':')
        print >> sys.stderr, "[%s] %s run initiated." % (time.ctime(), prog_name)
    else:
        sys.exit(("\nUsage: python %s <in.vcf> <inh.flt> <inh.type> <pat_id:flt> <fat_id:flt> <mot_id:flt>\n" % prog_name))
    # fi

    if pat_id == "-":
        sys.exit("ERROR: Patiend ID must be present.")

    pat_abgt_list = get_abgt_list(pat_flt)
    fat_abgt_list = get_abgt_list(fat_flt)
    mot_abgt_list = get_abgt_list(mot_flt)

    print_filtered_vcf(in_vcf, inh_type, pat_id, fat_id, mot_id,
                       pat_abgt_list, fat_abgt_list, mot_abgt_list)


