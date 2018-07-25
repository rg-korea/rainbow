############
##  MAIN  ##
############
from __future__ import print_function
import sys, re
from collections import defaultdict

prog_name = sys.argv[0].split('/')[-1]
if len(sys.argv) == 4:
    in_csv = sys.argv[1]
    in_omim = sys.argv[2] # genemap2.txt
    skip_if_no_omim_annot = int(sys.argv[3])
else:
    sys.exit("\nUsage: python %s <in.csv> <in.omim> <skip.if.no.omim.annotation[0/1]>\n" % prog_name)
# fi

class Phenotype:
    def __init__(self, info):
        # Leukodystrophy, hypomyelinating, 10, 616420 (3), Autosomal recessive
        self.info = info
        self.fields = self.info.split(', ')
        self.skip = False
        result = re.search("(.+), (\d+) \([1-4]\)", info)
        if not result:
            self.skip = True
        else:
            self.omim_desc, self.omim_id = result.groups()
            self.omim_desc = self.omim_desc.replace(", ", '; ')
            self.omim_link = "https://www.omim.org/entry/" + self.omim_id

def process_omim_genemap(in_omim):
    gene_to_phenotypes = defaultdict(list)
    for line in open(in_omim, "r"):
        if line.startswith("#"):
            continue
        field = line[:-1].split('\t')
        hgnc_symbol = field[8]
        if not hgnc_symbol: # empty column
            continue
        phenotypes = field[12].split('; ')
        if not phenotypes: # empty column
            continue
        flag_skip = True
        for info in phenotypes:
            if (not info or
                info.count('{') or
                info.count('[') or
                info.startswith('?')):
                continue
            else:
                phenotype = Phenotype(info)
                if phenotype.skip:
                    continue
                flag_skip = False
                gene_to_phenotypes[hgnc_symbol].append(phenotype)
    return gene_to_phenotypes

gene_to_phenotypes = process_omim_genemap(in_omim)

flag_header = True
for line in open(in_csv, "r"):
    if flag_header:
        flag_header = False
        delimitors = ',' * line.count(',')
        print(line.rstrip() + ',' + 'OMIM phenotype' + ',' + 'OMIM link')
        continue
    field = line.rstrip().split(',')
    hgnc_symbols = list()
    genes = field[5].split(';')
    for gene in genes:
        if gene not in hgnc_symbols:
            hgnc_symbols.append(gene)
    first_print = True
    phenotype_found = False
    for hgnc_symbol in hgnc_symbols:
        if hgnc_symbol in gene_to_phenotypes:
            phenotypes = gene_to_phenotypes[hgnc_symbol]
            phenotype_found = len(phenotypes)
            for phenotype in phenotypes:
                if first_print:
                    first_print = False
                    print(line.rstrip() + ',' + phenotype.omim_desc + ',' + phenotype.omim_link)
                else:
                    print(delimitors + ',' + phenotype.omim_desc + ',' + phenotype.omim_link)
    if not phenotype_found:
        if skip_if_no_omim_annot:
            continue
        else:
            print(line.rstrip())
