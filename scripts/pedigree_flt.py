## GLOBALS
from __future__ import print_function
DEBUG = True

## CUTOFFS
alt_cnt_cutoff = 5 #3 #5 #2
depth_cutoff = 10 #15 #10
alt_freq_cutoff = 0.15
chrom_ids = [str(x) for x in range(1,23)] + ['X', 'Y', 'M', "MT"]

###############
##  MODULES  ##
###############
def get_abgt_list(indiv_flt):
    if indiv_flt == "not_ref":
        return {"het1", "het2", "hom"}
    elif indiv_flt == "all":
        return {"ref", "het1", "het2", "hom"}
    elif indiv_flt == "het":
        return {"het1", "het2"}
    elif indiv_flt == "not_hom":
        return {"ref", "het1", "het2"}
    elif indiv_flt == "het1_hom":
        return {"het1", "hom"}
    elif indiv_flt == "ref_hom":
        return {"het1", "hom"}
    else:
        return {indiv_flt}
    # fi
# fed

def chrom_fit_pedigree(inh_type, chrom, chrom_ids):
    chrom_id = chrom.replace("chr",'')
    if chrom_id in chrom_ids:
        pass
    else:
        return False # remove contigs

    if inh_type.startswith('X'): # X-linked dominant/recessive
        if chrom_id.count('X'):
            return True
        else:
            return False
    elif inh_type.startswith('Y'): # Y-linked
        if chrom_id.count('Y'):
            return True
        else:
            return False
    elif inh_type.startswith('A'): # Autosomal dominant/homozygous/comp_hetero
        if (chrom_id.count('X') or chrom_id.count('Y') 
            or chrom_id.count('M') or chrom_id.count("MT")):
            return False # if is X, Y, or M/MT(mitochondrial)
        else:
            return True # autosomes
    elif inh_type == "MT": # MiTocondrial
        if chrom_id.count('M') or chrom_id.count("MT"):
            return True
        else:
            return False

    elif inh_type == "DN":
        return True
        
    else: # ERROR
        sys.exit("ERROR: inh_type: %s, chromosome:%s" % (inh_type, chrom))
    # fi
# fed

def proc_and_print_header(vcf_line, flt):
    field = vcf_line.strip().split('\t')
    flt.pat_idx = field.index(flt.pat_id)
    flt.fat_idx = field.index(flt.fat_id) if flt.fat_id != "-" else None
    flt.mot_idx = field.index(flt.mot_id) if flt.mot_id != "-" else None
    flt.have_both_parent_samples = (isinstance(flt.fat_idx, int) and 
                                    isinstance(flt.mot_idx, int))
    if flt.gene_count_phase:
        pass
    else:
        print(vcf_line.strip())

def get_alt_index(pat_genotype):
    if 0 in pat_genotype.ais:
        alt_idxs = list(filter(lambda x: x!=0, pat_genotype.ais))
    elif len(set(pat_genotype.ais)) == 1:
        alt_idxs = [pat_genotype.ais[0]]
    else:
        alt_idxs = pat_genotype.ais
    return alt_idxs

class Line:
    def __init__(self, vcf_line, flt):
        self.vcf_line = vcf_line
        self.field = vcf_line.rstrip().split('\t')
        self.chrom = self.field[0] # chr3
        self.pos = self.field[1] # 42251577 (str)
        self.chr_pos = "%s:%s" % (self.chrom, self.pos)
        self.ref = self.field[3]
        self.alts = self.field[4].split(',')
        self.alleles = [self.ref] + self.alts # 0:A[ref], 1:C[alt], 2:CT[alt], ...
        self.filter = set(self.field[6].split(';'))
        self.info = self.field[7]
        self.genotype_field = self.field[8] # GT:AD:DP:GQ:PL
        self.pat_genotype_info = self.field[flt.pat_idx]
        self.fat_genotype_info = None
        self.mot_genotype_info = None
        if flt.fat_idx:
            self.fat_genotype_info = self.field[flt.fat_idx]
        if flt.mot_idx:
            self.mot_genotype_info = self.field[flt.mot_idx]
        self.pflag = True # init

    def filter_patient_genotypes(self, flt):
        self.pat_genotype = Genotype(self, self.pat_genotype_info)
        self.fat_genotype = Genotype(self, self.fat_genotype_info)
        self.mot_genotype = Genotype(self, self.mot_genotype_info)

        if not self.pat_genotype.has_genotype: # has genotype?
            self.pflag = ('no genotype [%s]: %s' % (self.chr_pos,
                self.pat_genotype.has_genotype))
            return False
        if self.pat_genotype.agt == "ref": # has mutation?
            self.pflag = ('ref agt [%s]: %s' % (self.chr_pos, self.pat_genotype.agt))
            return False
        else:
            candidate_genes = list(filter(lambda x: x != '.',
                                   self.pat_genotype.candidate_alleles))
            self.pat_genotype.alt_idxs = get_alt_index(self.pat_genotype)
        if not chrom_fit_pedigree(flt.inh_type, self.chrom, flt.chrom_ids): # good chr?
            self.pflag = ('chrom not fit [%s]: %s' % (self.chr_pos, self.chrom))
            return False

        # good quality variant?
        if not (self.filter & {'PASS', '.'}): # good FILTER?
            self.pflag = ('FILTER bad [%s]: %s' % (self.chr_pos, self.filter))
            return False
        if ((self.pat_genotype.depth < flt.depth_cutoff) or
            (self.fat_genotype.has_genotype and self.fat_genotype.depth < flt.depth_cutoff) or
            (self.mot_genotype.has_genotype and self.mot_genotype.depth < flt.depth_cutoff)):
            self.pflag = ('depth low [%s]: pat %s / fat&mot not shown' % (self.chr_pos, self.pat_genotype.depth))
            return False
        else:
            self.pat_genotype.alt_cnt = max([self.pat_genotype.allele_depths[x] 
                                             for x in self.pat_genotype.alt_idxs])
            self.pat_genotype.alt_fraction = (float(self.pat_genotype.alt_cnt) 
                                              / self.pat_genotype.depth)
        if self.pat_genotype.alt_cnt < flt.alt_cnt_cutoff:
            self.pflag = ('alt_cnt low [%s]: %s' % (self.chr_pos, self.pat_genotype.alt_cnt))
            return False
        if self.pat_genotype.alt_fraction < flt.alt_freq_cutoff:
            self.pflag = ('alt_fraction low [%s]: %s' % (self.chr_pos, 
                self.pat_genotype.alt_fraction))
            return False

        # has functional candidate variant?
        if not candidate_genes:
            self.pflag = ('candidate_genes empty [%s]: %s' % (self.chr_pos, candidate_genes))
            return False

        # has fit pedigree?
        if self.pat_genotype.agt not in flt.pat_abgt_list:
            self.pflag = ('agt(%s) not in %s [%s]' % (
                self.pat_genotype.agt, flt.pat_abgt_list, self.chr_pos))
            return False
        if self.fat_genotype.agt not in flt.fat_abgt_list:
            self.pflag = ('agt(%s) not in %s [%s]' % (
                self.fat_genotype.agt, flt.fat_abgt_list, self.chr_pos))
            return False
        if self.mot_genotype.agt not in flt.mot_abgt_list:
            self.pflag = ('agt(%s) not in %s [%s]' % (
                self.mot_genotype.agt, flt.mot_abgt_list, self.chr_pos))
            return False

        # good inheritance status?
        if flt.inh_flt == 0: # inherited
            if self.pat_genotype.inheritance_status not in {'0', '1'}: 
                self.pflag = ('inheritance_status [%s] not in {0, 1}'
                    % self.chr_pos)
                return False

        # good compound hetero candidate?
        elif flt.inh_flt == 1: # comp het flag
            parent_inheritance_status = str(int(flt.have_both_parent_samples))
            if self.pat_genotype.inheritance_status not in {parent_inheritance_status, '3'}:
                self.pflag = ('inheritance_status [%s] %s not in {%s, 3}'
                    % (self.chr_pos, self.pat_genotype.inheritance_status, parent_inheritance_status))
                return False

            if flt.gene_count_phase: # gene count phase
                self.pat_genotype.get_parentage(self)
                for hgnc in candidate_genes: # add genes
                    for pat_ai in self.pat_genotype.ais:
                        if not pat_ai.is_variant: 
                            self.pflag = ('pat_ai.is_variant [%s] False'
                                % self.chr_pos)
                            continue
                        # add variant gene cnt
                        flt.var_gene_cnt[hgnc] += 1
                        # add var parentages to dict[hgnc]=set
                        flt.var_parentages_per_gene.setdefault(
                            hgnc, set()).add(pat_ai.parentage)
            else: # print phase
                len_comp_het_genes = [len(flt.var_parentages_per_gene[x])
                    for x in candidate_genes
                    if x in flt.var_parentages_per_gene]
                if flt.have_both_parent_samples: # from father & mother?
                    if ((not len_comp_het_genes) or # no variant allele parentage 
                        (max(len_comp_het_genes) < 2) or # parentage set size small
                        (min([flt.var_gene_cnt[x] 
                            for x in candidate_genes]) > 5)): # gene has too many variants
                        self.pflag = ('len_comp_het_genes = %s [%s] | flt.var_gene_cnt[list]=%s'
                            % (len_comp_het_genes, self.chr_pos, [flt.var_gene_cnt[x] for x in candidate_genes]))
                        return False
                else:
                    if ((not len_comp_het_genes) or # no variant allele parentage 
                        (max([flt.var_gene_cnt[x] 
                            for x in candidate_genes]) < 2) or # gene has too few variants
                        (min([flt.var_gene_cnt[x] 
                            for x in candidate_genes]) > 5)): # gene has too many variants
                        self.pflag = ('len_comp_het_genes = %s [%s] | flt.var_gene_cnt[list]=%s'
                            % (len_comp_het_genes, self.chr_pos, [flt.var_gene_cnt[x] for x in candidate_genes]))
                        return False
                        
        elif flt.inh_flt == 2: # de novo
            if self.pat_genotype.inheritance_status not in {'2', '3'}:
                self.pflag = ('inheritance_status [%s] not in {2, 3}'
                    % self.pat_genotype.inheritance_status)
                return False
        else:
            sys.exit("ERROR: Wrong inheritance status.")

        return True


class AlleleIndex(int):
    def __init__(self, value):
        self._value = value
        self.parentage = 00

class Genotype:
    def __init__(self, line, indiv_genotype_info):
        self.agt = '-' # init
        self.inheritance_status = '-' # init
        self.candidate_alleles = ['.']
        genotype_fields = line.genotype_field.split(':') # GT:AD:DP:GQ:PL -> GT, AD, DP, GQ, PL
        if (indiv_genotype_info and indiv_genotype_info != './.'): # if genotype info present
            self.indiv_genotype_infos = indiv_genotype_info.split(':')
            zipped_infos = zip(genotype_fields, self.indiv_genotype_infos)
            self.indiv_genotype_infos = {x[0]:x[1] for x in zipped_infos}
            self.gt = self.indiv_genotype_infos['GT'] # 1/2
            if self.gt == './.':
                self.has_genotype = False
            elif ('AD' in self.indiv_genotype_infos 
                  and self.indiv_genotype_infos['AD'].count('.')):
                self.has_genotype = False
            elif len(set(self.indiv_genotype_infos.keys()) 
                  & {'GT', 'AD', 'DP'}) != 3:
                self.has_genotype = False
            else:
                self.has_genotype = True
                self.ais = [AlleleIndex(x) for x in self.gt.split('/')]
                self.ad = self.indiv_genotype_infos['AD'] # 4,31,36
                self.allele_depths = [int(x) for x in self.ad.split(',')]
                self.dp = self.indiv_genotype_infos['DP'] # 71
                self.depth = int(self.dp)
                self.gq = self.indiv_genotype_infos['GQ'] # 99
                self.pl = self.indiv_genotype_infos['PL'] # 3203,1522,1509,1091,0,1213
                self.agt = self.indiv_genotype_infos['AG'] # ref, het1, het2, hom
                # 0=inherited, 1=comp het, 2=de novo, 3=comp het & de novo
                self.inheritance_status = self.indiv_genotype_infos['IH']
                # TP53=1 | TP53=1,TP53-1=1 | ...
                self.candidate_alleles = self.indiv_genotype_infos['GS'].split(',') 
        else:
            self.has_genotype = False

    def get_parentage(self, line):
        if self.has_genotype: # if pat data available
            for pat_ai in line.pat_genotype.ais: # init
                if pat_ai == 0: # if allele == ref
                    pat_ai.is_variant = False
                else: # if allele == var
                    pat_ai.is_variant = True
            if line.fat_genotype.has_genotype: # if fat data available, label parentage from fat
                fat_ai_set = set(line.fat_genotype.ais)
                for pat_ai in line.pat_genotype.ais:
                    if pat_ai in fat_ai_set:
                        pat_ai.parentage = 10 # found in fat
            if line.mot_genotype.has_genotype: # if mot data available, label parentage from mot
                mot_ai_set = set(line.mot_genotype.ais)
                for pat_ai in line.pat_genotype.ais:
                    if pat_ai in mot_ai_set:
                        pat_ai.parentage += 1 # found in mot
        else: # if pat data inavailable, do nothing
            pass

class FilterAttributes:
    def __init__(self, inh_flt, inh_type, pat_id, fat_id, mot_id, 
                 pat_flt, fat_flt, mot_flt, chrom_ids,
                 alt_cnt_cutoff, depth_cutoff, alt_freq_cutoff):
        self.inh_flt = int(inh_flt)
        assert self.inh_flt in {0, 1, 2}, "ERROR: inh_flt=%s" % self.inh_flt
        self.inh_type = inh_type
        self.pat_id = pat_id
        self.fat_id = fat_id
        self.mot_id = mot_id
        self.pat_abgt_list = get_abgt_list(pat_flt)
        self.fat_abgt_list = get_abgt_list(fat_flt)
        self.mot_abgt_list = get_abgt_list(mot_flt)
        self.chrom_ids = chrom_ids
        self.alt_cnt_cutoff = alt_cnt_cutoff
        self.depth_cutoff = depth_cutoff
        self.alt_freq_cutoff = alt_freq_cutoff
        self.var_parentages_per_gene = {}
        self.var_gene_cnt = defaultdict(int)
        self.gene_count_phase = True

############
##  MAIN  ##
############
    # for line end

if __name__ == "__main__":
    import sys, time, gzip
    from collections import defaultdict
    prog_name = sys.argv[0].split('/')[-1]
    if len(sys.argv) == 7:
        in_vcf = sys.argv[1]
        inh_flt = sys.argv[2]
        inh_type = sys.argv[3]
        pat_id, pat_flt = sys.argv[4].split(':')
        assert pat_id != "-", "ERROR: Patient ID must be present."
        fat_id, fat_flt = sys.argv[5].split(':')
        mot_id, mot_flt = sys.argv[6].split(':')
        flt = FilterAttributes(inh_flt, inh_type, pat_id, fat_id, mot_id,
                               pat_flt, fat_flt, mot_flt, chrom_ids,
                               alt_cnt_cutoff, depth_cutoff, alt_freq_cutoff)
        print("[%s] %s run initiated." % (time.ctime(), prog_name), file=sys.stderr)
    else:
        sys.exit(("\nUsage: python %s <in.vcf> <inh.flt> <inh.type> " % prog_name
            + "<pat_id:flt> <fat_id:flt> <mot_id:flt>\n"))
    # fi
  
    # Read phase, gene count phase
    if in_vcf.endswith('.gz'):
        vcf = gzip.open(in_vcf, "r")
    else:
        vcf = open(in_vcf, "r")
    for vcf_line in vcf:
        if vcf_line.startswith('##'):
            continue
        elif vcf_line.startswith('#'): # get sample genotype column index
            proc_and_print_header(vcf_line, flt)
            continue
        line = Line(vcf_line, flt)
        _ = line.filter_patient_genotypes(flt)
    vcf.close()

    # Print phase
    flt.gene_count_phase = False # print phase
    if in_vcf.endswith('.gz'):
        vcf = gzip.open(in_vcf, "r")
    else:
        vcf = open(in_vcf, "r")
    for vcf_line in vcf:
        if vcf_line.startswith('##'):
            print(vcf_line.rstrip())
            continue
        elif vcf_line.startswith('#'): # get sample genotype column index
            proc_and_print_header(vcf_line, flt)
            continue
        line = Line(vcf_line, flt)
        print_flag = line.filter_patient_genotypes(flt)
        if DEBUG:
            print(line.pflag, file=sys.stderr)
        if print_flag:
            print(line.vcf_line.rstrip())
    vcf.close()
