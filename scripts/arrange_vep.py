# Last update: July 1st 2016
# Author: Seongmin Choi <seongmin.choi@raregenomics.org>

############
##  INIT  ##
############

import sys, re, time

prog_name = sys.argv[0].split('/')[-1]
if len(sys.argv) == 4:
    in_vep = sys.argv[1]
    in_vcf = sys.argv[2]
    pat_id = sys.argv[3]
    print >> sys.stderr, "[%s] %s run initiated." % (time.ctime(), prog_name)
else:
    sys.exit("\nUsage: python %s <in.vep> <in.vcf> <pat_id>\n" % prog_name)
# fi

idx_s = 9 # column number from where SMP_ID starts


##############
##  MODULE  ##
##############

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

def get_vcf_data(in_vcf, pat_id, idx_s):
    var_info = dict()
    for line in open(in_vcf):
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
            
            formats = [ x.split(':') for x in field[idx_s:] ] # genotype information of chr:pos, for all samples
            smp_genos = { x: formats[smp_ids.index(x)] for x in smp_ids } # smp_genos[sample_id] = allele_genotype_information
            smp_genos["-"] = "0/1/2/3/4/5/6:0.0,0.0,0.0:0:0:0:0".split(':') # if sample not present
            # 0/1 -209.52,0.0,-196.62 1 99 152 76 het1 inherited
            depth = smp_genos[pat_id] [4]
            var_depth = smp_genos[pat_id] [5]
    
            if depth.count(',') and not alt_base.count(','):
                depth = depth.split(',')[-1]
                var_depth = var_depth.split(',')[-1]
            else:
                pass
            # fi
            var_tuple = ( ref_base, alt_base, smp_genos, depth, var_depth )
            var_info[chrpos] = var_tuple
        # fi
    # for line end
    
    return var_info
# fed

def print_output(var_info, in_vep): 
    printed_gene = dict()
    header = ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % 
              ("Genomic location", "Ref. allele", "Var. allele", "Gene symbol", "Consequence",
               "Amino acids", "Codons", "Read depth", "Var. allele count", "SIFT", "PolyPhen"))
    print header
    
    for line in open(in_vep):
        if line.startswith('#'): # skip comments
            continue
        
        field = line.strip().split('\t')
        if len(field) < 14: # skip improper lines
            continue
    
        (upvar, loc, al, gene_no, nm_id, feat, result, cdnap, cdsp, 
         protp, aas, codons, known_var, extra) = field
    
        if not nm_id.startswith("NM_"): # skip non NM's
            continue
    
        chrpos = get_chrpos(upvar)
        #ref_al, var_al = upvar.split('_')[-1].split('/')
        chrome, pos = chrpos.split(':')
        var_tuple = var_info[chrpos]
        ( ref_base, alt_base, smp_genos, depth, var_depth ) = var_tuple
    
        # STRAND=-1;SYMBOL=TBC1D32;SIFT=tolerated(0.2);PolyPhen=benign(0.007)
        gene_sym = re.search("SYMBOL=([A-Z\d]+)", extra).groups()[0]
        if gene_sym in printed_gene: # skip if geneSym is already printed
            if printed_gene[gene_sym] == chrpos:
                continue
            else:
                printed_gene[gene_sym] = chrpos
        else:
            printed_gene[gene_sym] = chrpos
        #if gene_sym=="CEP350": ##@##
        #    print line
        sift_re = re.search("SIFT=(.+?);", extra)
        polyphen_re = re.search("PolyPhen=(.+?)\)", extra)
        sift = sift_re.groups()[0] if sift_re else ""
        polyphen = polyphen_re.groups()[0]+')' if polyphen_re else ""
        
        print_line = ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                      (chrpos, ref_base, alt_base, gene_sym, result, 
                       aas, codons, depth, var_depth, sift, polyphen))
        print print_line
    
    # for line end
# fed


############
##  MAIN  ##
############

var_info = get_vcf_data(in_vcf, pat_id, idx_s)
print_output(var_info, in_vep)
   
