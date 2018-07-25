############
##  INIT  ##
############
from __future__ import print_function
idx_s = 9 # column number from where SMP_ID starts
delim = ',' # deliminator of output result

##############
##  MODULE  ##
##############
class Annotation:
    def __init__(self, annotation_string):
        # 0 1                2        3    4    5          6           7              8     9         10          11        12       13       14 15
        # A|missense_variant|MODERATE|KAZN|KAZN|transcript|NM_201628.2|protein_coding|12/15|c.1796G>A|p.Arg599Gln|2077/6022|1796/2328|599/775||
        separated_annotation = annotation_string.split('|')
        self.allele = separated_annotation[0]
        self.effect = separated_annotation[1]
        self.effects = set(self.effect.split('&'))
        self.impact = separated_annotation[2]
        self.hgnc_symbol = separated_annotation[3]
        self.refseq_id = separated_annotation[6]
        self.exon = separated_annotation[8]
        self.cdna_var = separated_annotation[9]
        self.aa_var = separated_annotation[10]
        self.cdna_pos = separated_annotation[11]
        self.cds_pos = separated_annotation[12]
        self.aa_pos = separated_annotation[13]
        if self.impact == 'HIGH' or self.impact == 'MODERATE':
            self.flag_save = True 
        else:
            self.flag_save = False

class AlleleIndex(int):
    def __init__(self, value):
        self._value = value
 
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
                self.var_ais = [AlleleIndex(x) for x in self.ais if x != 0]
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

class VcfAttributes:
    def __init__(self, pat_id):
        self.pat_id = pat_id
 
class Line:
    def __init__(self, vcf_line, vcf_attr):
        self.vcf_line = vcf_line
        self.field = vcf_line.rstrip().split('\t')
        self.chrom = self.field[0] # chr3
        self.pos = self.field[1] # 42251577 (str)
        self.chr_pos = "%s:%s" % (self.chrom, self.pos)
        self.ref = self.field[3]
        self.alt = self.field[4]
        self.alts = self.alt.split(',')
        self.alleles = [self.ref] + self.alts # 0:A[ref], 1:C[alt], 2:CT[alt], ...
        self.filter = set(self.field[6].split(';'))
        self.info = self.field[7]
        self.genotype_field = self.field[8] # GT:AD:DP:GQ:PL
        self.pat_genotype_info = self.field[vcf_attr.pat_idx]
        self.pflag = True # init
        self.pat_genotype = Genotype(self, self.pat_genotype_info)

    def process_snpeff_annotation_in_info_field(self):
        assert 'ANN=' in self.info, "ERROR: %s does not have ANN" % infos
        infos = self.info.split(';')
        annotation_info = [x[4:] for x in infos if x.startswith('ANN=')][0]
        annotation_strings = annotation_info.split(',')
        self.annotation_by_allele = defaultdict(list)
        self.candidate_present = False # default false for candidate presence
        self.var_alleles = {self.alleles[x] for x in self.pat_genotype.var_ais}
        self.var_ads = {self.pat_genotype.allele_depths[x] 
            for x in self.pat_genotype.var_ais}
        self.ad_sum = sum(self.pat_genotype.allele_depths)
        self.var_vafs = {self.pat_genotype.allele_depths[i]/float(self.ad_sum)
            for i in self.pat_genotype.var_ais}
        self.effects = []
        self.exons = []
        self.cdna_vars = []
        self.aa_vars = []
        self.cdna_poss = []
        self.cds_poss = []
        self.aa_poss = []
        self.hgnc_symbols = []
        self.refseq_ids = []
        for annotation_string in annotation_strings:
            annotation = Annotation(annotation_string)
            if (annotation.flag_save and
                (annotation.allele in self.var_alleles) and
                annotation.refseq_id.startswith("NM_")):
                self.candidate_present = True
                self.effects.append(annotation.effect)
                self.exons.append(annotation.exon)
                self.cdna_vars.append(annotation.cdna_var)
                self.aa_vars.append(annotation.aa_var)
                self.cdna_poss.append(annotation.cdna_pos)
                self.cds_poss.append(annotation.cds_pos)
                self.aa_poss.append(annotation.aa_pos)
                self.hgnc_symbols.append(annotation.hgnc_symbol)
                self.refseq_ids.append(annotation.refseq_id)
                self.annotation_by_allele[annotation.allele].append(annotation)
            else:
                continue

def proc_and_print_header(vcf_line, vcf_attr):
    field = vcf_line.strip().split('\t')
    vcf_attr.pat_idx = field.index(vcf_attr.pat_id)


############
##  MAIN  ##
############
if __name__ == "__main__":
    import sys, time
    from collections import defaultdict
    prog_name = sys.argv[0].split('/')[-1]
    if len(sys.argv) == 3:
        in_vcf = sys.argv[1]
        pat_id = sys.argv[2]
        vcf_attr = VcfAttributes(pat_id)
        print("[%s] %s run initiated." % (time.ctime(), prog_name), file=sys.stderr)
    else:
        sys.exit("\nUsage: python %s <in.vcf> <pat_id>\n" % prog_name)

    headers = ("Genomic coordinate", "Reference", "Variants", "VAFs", "Depth",
               "HGNC Symbols", "RefSeq IDs", "Effects", "Exon numbers", "cDNA annotations", "AA annotations",
               "cDNA positions", "CDS positions", "AA positions")
    print(delim.join(headers))

    for vcf_line in open(in_vcf, "r"):
        if vcf_line.startswith('##'):
            continue
        elif vcf_line.startswith('#'): # get sample genotype column index
            proc_and_print_header(vcf_line, vcf_attr)
            continue
        line = Line(vcf_line, vcf_attr)
        line.process_snpeff_annotation_in_info_field()
        genomic_coord = line.chrom + ':' + line.pos
        var_alleles = ';'.join(line.var_alleles)
        hgnc_symbols = ';'.join(line.hgnc_symbols)
        refseq_ids = ';'.join(line.refseq_ids)
        effects = ';'.join(line.effects)
        exons = ';'.join(line.exons)
        cdna_vars = ';'.join(line.cdna_vars)
        aa_vars = ';'.join(line.aa_vars)
        cdna_poss = ';'.join(line.cdna_poss)
        cds_poss = ';'.join(line.cds_poss)
        aa_poss = ';'.join(line.aa_poss)
        depth = str(line.ad_sum)
        vafs = ';'.join(["%.2f"%x for x in line.var_vafs])
        outputs = (genomic_coord, line.ref, var_alleles, vafs, depth,
            hgnc_symbols, refseq_ids, effects, exons, cdna_vars, aa_vars, 
            cdna_poss, cds_poss, aa_poss)
        assert len(headers) == len(outputs), "ERROR: header len != output len"
        fmt = (("{}"+delim) * len(outputs)).rstrip(delim)
        print(fmt.format(*outputs))

