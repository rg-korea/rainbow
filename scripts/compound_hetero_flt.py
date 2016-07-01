# Last update: July 1st 2016
# Author: Seongmin Choi <seongmin.choi@raregenomics.org>

############
##  INIT  ##
############

import os, sys, re

if len(sys.argv) == 6:
    in_vcf = sys.argv[1]
    in_vep = sys.argv[2]
    pat_id = sys.argv[3]
    fat_id = sys.argv[4]
    fat_id = fat_id if fat_id!="-" else "-f"
    mot_id = sys.argv[5]
    mot_id = mot_id if mot_id!="-" else "-m"
else:
    prog_name = sys.argv[0].split('/')[-1]
    sys.exit("\nUsage: python2.7 %s <in.trio.vcf> <in.vep> <pat_id> <fat_id> <mot_id>\n" % prog_name)
# fi

idx_s = 9 # where sample names start, from the final comment line
mutated = ["het1", "het2", "hom"]
high_impact = ["transcript_ablation","splice_acceptor_variant","splice_donor_variant",
    "stop_gained","frameshift_variant","stop_lost","start_lost","transcript_amplification",
    "inframe_insertion","inframe_deletion","missense_variant","protein_altering_variant",
    "incomplete_terminal_codon_variant","synonymous_variant","coding_sequence_variant"]


###############
##  MODULES  ##
###############

# Convert VEP position to VCF position
def get_chrpos(upvar):
    chr_no, vep_pos, ref_alt = upvar.split('_')
    chrom = "chr%s" % chr_no
    ref = ref_alt.split('/')[0]
    alts = ref_alt.split('/')
    
    if ref == '-' or (len(ref) < max( [len(x) for x in alts] )):
        one_pos = int(vep_pos) - 1
    elif '-' in alts or (len(ref) > min( [len(x) for x in alts] )):
        one_pos = int(vep_pos) - 1
    else:
        one_pos = int(vep_pos)
    
    chrpos = "%s:%s" % (chrom, one_pos)
    return chrpos
# fed 

# Proc VEP
def get_vep_data(in_vep):
    vep_data = dict()
    for line in open(in_vep, "r"):
        if line.startswith('#'): # skip comments
            continue
        
        field = line.strip().split('\t')
        if len(field) < 14: # skip improper lines
            continue
    
        (upvar, loc, al, gene_no, nm_id, feat, result, cdnap, cdsp, 
         protp, aas, codons, known_var, extra) = field

        chrpos = get_chrpos(upvar)    
    
        if not nm_id.startswith("NM_"): # skip non NM's
            continue
    
        gene_sym = re.search("SYMBOL=([A-Z\d]+)", extra).groups()[0]
        if chrpos not in vep_data:
            vep_data[chrpos] = {nm_id: (gene_sym, set(result.split(',')))}
        else:
            vep_data[chrpos][nm_id] = (gene_sym, set(result.split(',')))
        # fi
    # for line end

    return vep_data
# fed
    
def is_comp_het(field, idx_s, smp_ids, 
                pat_id, fat_id, mot_id, mutated):
    # 1/0 : -6.79,0.0,-37.39 : 11 : 68 : 19 : 6 : het1 : inherited
    hets = [ x.split(':')[6] for x in field[idx_s:] ] # heterogeneity type of chr:pos, for all samples
    pat_het = hets[smp_ids.index(pat_id)]
    fat_het = hets[smp_ids.index(fat_id)] if fat_id != "-f" else False
    mot_het = hets[smp_ids.index(mot_id)] if mot_id != "-m" else False

    is_from_fat = (pat_het == "het1" and fat_het in mutated and mot_het == "ref")
    is_from_mot = (pat_het == "het1" and mot_het in mutated and fat_het == "ref")

    if fat_id == "-f" or mot_id == "-m":
        return (pat_het == "het1"), "-f", "-m"
    else:
        return (is_from_fat ^ is_from_mot), is_from_fat, is_from_mot # if comp het candidate
# fed

def add_up_var_data( nm_id, var_info, info_tuple, 
                     pat_id, fat_id, mot_id, is_from_fat, is_from_mot ):
    if nm_id not in var_info:
        var_info[nm_id] = dict()
    
    if pat_id in var_info[nm_id]:
        var_info[nm_id][pat_id].append( info_tuple ) # (comp_het, chrpos, gene_sym, results)
    else:
        var_info[nm_id][pat_id] = [ info_tuple ]
    # fi
    
    if is_from_fat == True: # if mut is from father
        if fat_id in var_info[nm_id]:
            var_info[nm_id][fat_id].append( info_tuple )
        else:
            var_info[nm_id][fat_id] = [ info_tuple ]
        # fi
    elif is_from_mot == True: # else if mut is from mother
        if mot_id in var_info[nm_id]:
            var_info[nm_id][mot_id].append( info_tuple )
        else:
            var_info[nm_id][mot_id] = [ info_tuple ]
        # fi
    else:
        return var_info
    # fi

    return var_info
# fed

# Proc VCF
def get_vcf_data(in_vcf, idx_s, pat_id, fat_id, mot_id, mutated, high_impact):
    var_info = dict()
    for line in open(in_vcf, "r"):
        if line.startswith("##"):
            continue
        elif line.startswith('#'):
            smp_ids = line.strip().split('\t') [idx_s:] # get sample names
            continue
        else:
            field = line.strip().split('\t')
            chrom = field[0] # chromosome
            if len(chrom) > 5: # remove chr19_gl000209_random etc.
                continue
            one_pos = int(field[1]) # one-based position
            chrpos = "%s:%s" % (chrom, one_pos) # chr:pos allele position

            comp_het, is_from_fat, is_from_mot = is_comp_het(field, idx_s, smp_ids,
                                                             pat_id, fat_id, mot_id, mutated)
            if chrpos in vep_data:
                nm_ids = vep_data[chrpos].keys() # get gene symbols
            else:
                continue
    
            for nm_id in nm_ids:
                gene_sym, results = vep_data[chrpos][nm_id] # get mutation consequence
                if not results.intersection(high_impact): # if consequence is not impactful
                    continue
                info_tuple = (comp_het, chrpos, gene_sym, results)

                var_info = add_up_var_data( nm_id, var_info, info_tuple, 
                    pat_id, fat_id, mot_id, is_from_fat, is_from_mot )
            # for nm_id end
        # fi
    # for line end

    return var_info
# fed

# Is VCF line compound heterozygotic?
def is_comp_het_line(vep_data, var_info, chrpos, pat_id, fat_id, mot_id):
    comp_het_line, gene_sym = False, False # init
    if chrpos in vep_data:
        nm_ids = vep_data[chrpos].keys() # get gene symbols
    else:
        #print >> sys.stderr, ("LOG: %s not in VEP" % chrpos)
        return False, False

    for nm_id in nm_ids:
        if nm_id in var_info: # if saved
            pat_muts = filter(lambda x:x[0], var_info[nm_id][pat_id])
            gene_sym = ",".join([x[2] for x in pat_muts])

            if fat_id in var_info[nm_id]:
                fat_muts = filter(lambda x:x[0], var_info[nm_id][fat_id])
            else:
                fat_muts = []
            # fi
            if mot_id in var_info[nm_id]:
                mot_muts = filter(lambda x:x[0], var_info[nm_id][mot_id])
            else:
                mot_muts = []
            # fi
        # fi
        else:
            continue

        if len(pat_muts) == 2 and len(fat_muts) == 1 and len(mot_muts) == 1:
            comp_het_line = True
        elif len(pat_muts) == 2 and fat_id == "-f" and mot_id == "-m":
            comp_het_line = True
        else:
            continue
        # fi
    # for nm_id end
    return comp_het_line, gene_sym
# fed

# Print compound hetiero variants
def print_comp_hets(in_vcf, vep_data, var_info, pat_id, fat_id, mot_id):
    for line in open(in_vcf, "r"): 
        if line.startswith('#'): # skip comments
            print line.strip()
            continue
    
        field = line.strip().split('\t')
        chrom = field[0] # chromosome
        if len(chrom) > 5: # remove chr19_gl000209_random etc.
            continue
        one_pos = int(field[1]) # one-based position
        chrpos = "%s:%s" % (chrom, one_pos) # chr:pos allele position
        
        comp_het_line, gene_sym = is_comp_het_line(vep_data, var_info, chrpos, pat_id, fat_id, mot_id) 

        if comp_het_line:
            field[7] += ";GS=%s" % gene_sym
            print '\t'.join(field)
        else:
            continue
        # fi
    # for line end


############
##  MAIN  ##
############
vep_data = get_vep_data(in_vep)
var_info = get_vcf_data(in_vcf, idx_s, pat_id, fat_id, mot_id, mutated, high_impact)
print_comp_hets(in_vcf, vep_data, var_info, pat_id, fat_id, mot_id)

