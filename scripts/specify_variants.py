###############
##  GLOBALS  ##
###############
from __future__ import print_function
idx_s = 9 # column number from where sample info (e.g. lca_IV-1) starts
max_alt_other_for_de_novo = 0 # max alt_cnt allowed for the patient_var_pos in parents
max_alt_freq_other_for_de_novo = 0.02 # max alt_freq allowed for the patient_var_pos 
                                      # in parents where the depth of parent_alt_pos 
                                      # is over min_depth_other_for_alt_de_novo
min_depth_other_for_alt_de_novo = 50 # cutoff for high-depth alt_pos in parents, 
                                     # to determine if the alt_freq in parents are low enough
min_depth_other_for_de_novo = 10 # min_depth for parent_alt_pos
de_novo_flts = (max_alt_other_for_de_novo, # tuple of filter cutoffs
                max_alt_freq_other_for_de_novo, 
                min_depth_other_for_de_novo,
                min_depth_other_for_alt_de_novo)

DEBUG = False

###############
##  MODULES  ##
###############
def proc_and_print_header(vcf_line):
    field = vcf_line.strip().split('\t')
    pat_ids = [x for x in field[idx_s:] if (x!=fat_id and x!=mot_id)]
    pat_idxs = [field.index(x) for x in pat_ids]
    fat_idx = field.index( fat_id ) if fat_id != "-" else None
    mot_idx = field.index( mot_id ) if mot_id != "-" else None
    candidate_cnt_per_gene = {}
    for indiv_idx in pat_idxs + [fat_idx, mot_idx]:
        candidate_cnt_per_gene[indiv_idx] = defaultdict(int) # init
    print('##FORMAT=<ID=AG,Number=.,Type=String,Description='
          '"Abstract genotype: -, ref, het1, het2, hom">')
    print('##FORMAT=<ID=IH,Number=.,Type=Int,Description="'
          'Inheritence status: 0=inherited, '
          '1=compound heterozygous candidate, '
          '2=de novo candidate, '
          '3=compound heterozygous candidate and de novo candidate">')
    print(vcf_line.strip())
    return pat_idxs, fat_idx, mot_idx, candidate_cnt_per_gene

class Annotation:
    def __init__(self, annotation_string):
        separated_annotation = annotation_string.split('|')
        self.allele = separated_annotation[0]
        self.effect = separated_annotation[1]
        self.effects = set(self.effect.split('&'))
        self.impact = separated_annotation[2]
        self.hgnc_symbol = separated_annotation[3]
        self.exon_intron_rank = separated_annotation[4]
        self.cdna_var = separated_annotation[9]
        self.aa_var = separated_annotation[10]
        self.cdna_pos = separated_annotation[11]
        self.cds_pos = separated_annotation[12]
        self.aa_pos = separated_annotation[13]
        if self.impact == 'HIGH' or self.impact == 'MODERATE':
            self.flag_save = True 
        else:
            self.flag_save = False

class Line:
    def __init__(self, vcf_line, pat_idxs, fat_idx, mot_idx):
        self.vcf_line = vcf_line
        self.pat_idxs = pat_idxs[:]
        self.fat_idx = fat_idx
        self.mot_idx = mot_idx
        self.field = vcf_line.rstrip().split('\t')
        self.chrom = self.field[0] # chr3
        self.pos = self.field[1] # 42251577 (str)
        self.ref = self.field[3]
        self.alts = self.field[4].split(',')
        self.ref_alt_alleles = [self.ref] + self.alts # 0:A[ref], 1:C[alt], 2:CT[alt], ...
        self.info = self.field[7]
        self.genotype_field = self.field[8] # GT:AD:DP:GQ:PL
        self.pat_genotype_infos = [self.field[pat_idx] for pat_idx in pat_idxs]
        self.fat_genotype_info = None
        self.mot_genotype_info = None
        if fat_idx:
            self.fat_genotype_info = self.field[fat_idx]
        if mot_idx:
            self.mot_genotype_info = self.field[mot_idx]

    def process_snpeff_annotation_in_info_field(self):
        assert 'ANN=' in self.info, "ERROR: %s does not have ANN" % infos
        infos = self.info.split(';')
        annotation_info = [x[4:] for x in infos if x.startswith('ANN=')][0]
        annotation_strings = annotation_info.split(',')
        self.annotation_by_allele = defaultdict(list)
        self.candidate_present = False # default false for candidate presence
        for annotation_string in annotation_strings:
            # A|missense_variant|MODERATE|KAZN|KAZN|transcript|NM_201628.2|protein_coding|12/15|c.1796G>A|p.Arg599Gln|2077/6022|1796/2328|599/775||
            annotation = Annotation(annotation_string)
            if annotation.flag_save:
                self.candidate_present = True
                self.annotation_by_allele[annotation.allele].append(annotation)
            else:
                continue

    def process_family_genotypes(self, candidate_cnt_per_gene):
        self.fat_genotype = Genotype(self, self.fat_genotype_info, self.fat_idx)
        self.mot_genotype = Genotype(self, self.mot_genotype_info, self.mot_idx)
        self.pat_genotypes = []
        for i, pat_genotype_info in enumerate(self.pat_genotype_infos):
            pat_genotype = Genotype(self, pat_genotype_info, self.pat_idxs[i])
            self.pat_genotypes.append(pat_genotype)
            if pat_genotype.has_genotype:
                pat_genotype.get_abstract_genotype()
                pat_genotype.get_inheritance_status(line, pat_genotype, de_novo_flts)
                pat_genotype.get_affected_gene(self, candidate_cnt_per_gene)
                if DEBUG: print("pat cand:", self.chrom, self.pos, pat_genotype.agt, 
                    pat_genotype.ais, pat_genotype.depth, pat_genotype.allele_depths, 
                    "inh_state:", pat_genotype.inheritance_status, 
                    pat_genotype.candidate_alleles) 
            else: 
                if DEBUG: print("pat has no genotype...", pat_genotype.has_genotype)
        if self.fat_genotype.has_genotype:
            self.fat_genotype.get_abstract_genotype()
            self.fat_genotype.get_affected_gene(self, candidate_cnt_per_gene)
            if DEBUG: print("fat_cand:", self.fat_genotype.candidate_alleles)
        if self.mot_genotype.has_genotype:
            self.mot_genotype.get_abstract_genotype()
            self.mot_genotype.get_affected_gene(self, candidate_cnt_per_gene)
            if DEBUG: print("mot_cand:", self.mot_genotype.candidate_alleles)

    def print_line_with_additional_attributes(self):
        new_field = self.field[:]
        new_field[8] = self.genotype_field + ":" + "AG:IH:GS" # abst gt, inh status, gene symbol
        for i, pat_idx in enumerate(self.pat_idxs):
            pat_genotype = self.pat_genotypes[i]
            if pat_genotype.has_genotype:
                gs_subfield_str = '.'
                if pat_genotype.candidate_alleles:
                    gs_subfield_str = ','.join(list(pat_genotype.candidate_alleles))
                new_field[pat_idx] = ':'.join([self.field[pat_idx],
                                               pat_genotype.agt,
                                               str(pat_genotype.inheritance_status),
                                               gs_subfield_str])
        if self.fat_genotype.has_genotype:
            gs_subfield_str = '.'
            if self.fat_genotype.candidate_alleles:
                gs_subfield_str = ','.join(list(self.fat_genotype.candidate_alleles))
            new_field[self.fat_idx] = ':'.join([self.field[self.fat_idx],
                                                self.fat_genotype.agt,
                                                str(self.fat_genotype.inheritance_status),
                                                gs_subfield_str])
        if self.mot_genotype.has_genotype:
            gs_subfield_str = '.'
            if self.mot_genotype.candidate_alleles:
                gs_subfield_str = ','.join(list(self.mot_genotype.candidate_alleles))
            new_field[self.mot_idx] = ':'.join([self.field[self.mot_idx],
                                                self.mot_genotype.agt,
                                                str(self.mot_genotype.inheritance_status),
                                                gs_subfield_str])
        new_vcf_line = '\t'.join(new_field)
        print(new_vcf_line)
            
class AlleleIndex(int):
    def __init__(self, value):
        self._value = value
        # 10: found in fat, 01: found in mot, 11: found in both, 00: found in none
        self.parentage = 00 

def is_pat_alt_not_in_indiv(pat_genotype, indiv_genotype):
    not_in_indiv = False
    for pat_ai in pat_genotype.ais: # 0, 1 or 1, 2, ...
        if pat_ai == 0: # ignore patient ref allele
            continue
        elif pat_ai in indiv_genotype.ais:
            return False
        else:
            not_in_indiv = True
        # fi
    # for end

def is_indiv_alt_cnt_infrequent(indiv_genotype, de_novo_flts):
    (max_alt_other_for_de_novo, 
     max_alt_freq_other_for_de_novo, 
     min_depth_other_for_de_novo,
     min_depth_other_for_alt_de_novo) = de_novo_flts

    if indiv_genotype.depth < min_depth_other_for_de_novo:
        return False # low depth: not deep enough to ensure de novo

    max_alt_cnt = max(indiv_genotype.allele_depths[1:]) # max var allele depth
    max_alt_freq = float(max_alt_cnt) / indiv_genotype.depth

    if max_alt_cnt > max_alt_other_for_de_novo: # indiv has even 1 alt read
        # low freq & sufficient depth?
        if (max_alt_freq < max_alt_freq_other_for_de_novo and
            indiv_genotype.depth >= min_depth_other_for_alt_de_novo): 
            return True
        else:
            return False # not enough evidence to ensure de novo
        # fi
    else:
        return True # indiv has no alt read at patient alt position

def init_genotype(indiv_genotype, line, indiv_genotype_info, indiv_idx):
    genotype_fields = line.genotype_field.split(':') # GT:AD:DP:GQ:PL -> GT, AD, DP, GQ, PL
    if (indiv_genotype_info and indiv_genotype_info != './.'): # if genotype info present
        indiv_genotype.indiv_genotype_infos = indiv_genotype_info.split(':')
        zipped_infos = zip(genotype_fields, indiv_genotype.indiv_genotype_infos)
        indiv_genotype.indiv_genotype_infos = {x[0]:x[1] for x in zipped_infos}
        indiv_genotype.gt = indiv_genotype.indiv_genotype_infos['GT'] # 1/2
        indiv_genotype.indiv_idx = indiv_idx # unique id per indiv
        if indiv_genotype.gt == './.':
            indiv_genotype.has_genotype = False
        elif ('AD' in indiv_genotype.indiv_genotype_infos 
              and indiv_genotype.indiv_genotype_infos['AD'].count('.')):
            indiv_genotype.has_genotype = False
        elif len(set(indiv_genotype.indiv_genotype_infos.keys()) 
              & {'GT', 'AD', 'DP'}) != 3:
            indiv_genotype.has_genotype = False
        else:
            indiv_genotype.has_genotype = True
            indiv_genotype.ais = [AlleleIndex(x) for x in indiv_genotype.gt.split('/')]
            try:
                indiv_genotype.ad = indiv_genotype.indiv_genotype_infos['AD'] # 4,31,36
            except:
                sys.exit("indiv_genotype.indiv_genotype_infos = %s" % indiv_genotype.indiv_genotype_infos)
            indiv_genotype.allele_depths = [int(x) for x in indiv_genotype.ad.split(',')]
            indiv_genotype.dp = indiv_genotype.indiv_genotype_infos['DP'] # 71
            indiv_genotype.depth = int(indiv_genotype.dp)
            indiv_genotype.gq = indiv_genotype.indiv_genotype_infos['GQ'] # 99
            indiv_genotype.pl = indiv_genotype.indiv_genotype_infos['PL'] # 3203,1522,1509,1091,0,1213
            indiv_genotype.agt = '.' # init
            indiv_genotype.inheritance_status = '.' # init
            indiv_genotype.candidate_alleles = set() # init
    else:
        indiv_genotype.has_genotype = False

class Genotype:
    def __init__(self, line, indiv_genotype_info, indiv_idx):
        init_genotype(self, line, indiv_genotype_info, indiv_idx)

    def get_abstract_genotype(self):
        f_name = inspect.currentframe().f_code.co_name
        assert self.has_genotype, ("ERROR[%s]: "
            "%s does not have genotype." % (f_name, indiv_genotype_info))
        ai_set = set(self.ais)
        if self.gt == '0/0':
            self.agt = 'ref'
        elif self.ais.count(0) == 1: # only one 0[ref] <- 0/1, 0/2, ...
            self.agt = 'het1'
        elif len(ai_set) == 1: # no 0 and set size 1 <- 1/1, 2/2, ...
            self.agt = 'hom'
        elif len(ai_set) == 2: # no 0 and set size 2 <- 1/2, 2/4, ...
            self.agt = 'het2'
        else:
            sys.exit("ERROR: %s not classified as abstract." % self.gt)

    def get_inheritance_status(self, line, pat_genotype, de_novo_flts):
        if self.has_genotype: # if pat data available
            # inheritance status -> {0: inherited,
            #                        1: compound het candidate
            #                        2: de novo candidate
            #                        3: comp het candidate & de novo
            self.inheritance_status = 0 # default: inherited
            fat_has_no_alt = False # init
            mot_has_no_alt = False # init
            for pat_ai in pat_genotype.ais: # init
                if pat_ai == 0: # if allele == ref
                    pat_ai.is_variant = False
                else: # if allele == var
                    pat_ai.is_variant = True
            if line.fat_genotype.has_genotype: # if fat data available, label parentage from fat
                fat_has_no_alt = is_indiv_alt_cnt_infrequent(line.fat_genotype, de_novo_flts)
                fat_ai_set = set(line.fat_genotype.ais)
                for pat_ai in pat_genotype.ais:
                    if pat_ai in fat_ai_set:
                        pat_ai.parentage = 10 # found in fat
            if line.mot_genotype.has_genotype: # if mot data available, label parentage from mot
                mot_has_no_alt = is_indiv_alt_cnt_infrequent(line.mot_genotype, de_novo_flts)
                mot_ai_set = set(line.mot_genotype.ais)
                for pat_ai in pat_genotype.ais:
                    if pat_ai in mot_ai_set:
                        pat_ai.parentage += 1 # found in fat
            # comp het: pat.agt=het1 & ref found in parent, var can be de novo
            if pat_genotype.agt == 'het1':
                ref_ai = [x for x in pat_genotype.ais if x == 0][0]
                var_ai = [x for x in pat_genotype.ais if x != 0][0]
                # if ref from fat, var from mot or is de novo, vice versa
                if ((ref_ai.parentage in {10, 1, 11}) and
                    (var_ai.parentage in {1, 10, 00}) and
                    ((ref_ai.parentage + var_ai.parentage) in
                    {11, 12, 21, 10, 1})):
                    self.inheritance_status += 1 # compound het candidate
            # de novo: one of var allele is not found in both & parents don't have var
            if 00 in {x.parentage for x in pat_genotype.ais if x.is_variant}:
                if fat_has_no_alt and mot_has_no_alt:
                    self.inheritance_status += 2 # de novo candidate
        else: # if pat data inavailable, do nothing
            pass

    def get_affected_gene(self, line, candidate_cnt_per_gene):
        var_allele_idxs = filter(lambda x: x!=0, self.ais)
        add_candidate_var_to_gene = False
        affected_genes_per_variant = set() # all gene symbol annot for this var
        for var_allele_idx in var_allele_idxs:
            var_allele = line.ref_alt_alleles[var_allele_idx] # e.g. ACC
            if var_allele in line.annotation_by_allele: # var allele in dict
                candidate_annotations = line.annotation_by_allele[var_allele]
                affected_genes = set([x.hgnc_symbol for x in candidate_annotations])
                affected_genes_per_variant.update(affected_genes)
                add_candidate_var_to_gene = True
            else:
                continue
        for hgnc in affected_genes_per_variant:
            candidate_cnt_per_gene[self.indiv_idx][hgnc] += 1
            self.candidate_alleles.add(hgnc) # save gene

############
##  MAIN  ##
############
if __name__ == "__main__":
    import sys, time, inspect
    from collections import defaultdict
    prog_name = sys.argv[0].split('/')[-1]
    if len(sys.argv) == 4:
        in_vcf = sys.argv[1]
        fat_id = sys.argv[2]
        mot_id = sys.argv[3]
        print("[%s] %s run initiated." % (time.ctime(), prog_name), file=sys.stderr)
    else:
        sys.exit(("\nUsage: python %s <in.vcf> <fat_id> <mot_id>\n" % prog_name))

    for vcf_line in open(in_vcf):
        if vcf_line.startswith("##"): # skip header
            print(vcf_line.strip())
            continue
        elif vcf_line.startswith('#'): # get sample genotype column index
            (pat_idxs, fat_idx, mot_idx, 
             candidate_cnt_per_gene) = proc_and_print_header(vcf_line)
            continue
        line = Line(vcf_line, pat_idxs, fat_idx, mot_idx)
        line.process_snpeff_annotation_in_info_field()
        line.process_family_genotypes(candidate_cnt_per_gene)
        line.print_line_with_additional_attributes()
    # for line end

