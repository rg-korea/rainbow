# Last update: July 1st 2016
# Author: Seongmin Choi <seongmin.choi@raregenomics.org>

############
##  INIT  ##
############

import sys, marshal, bitarray, time

prog_name = sys.argv[0].split('/')[-1]
if len(sys.argv) == 3:
    in_vcf = sys.argv[1]
    in_bad_db = sys.argv[2]
    print >> sys.stderr, "[%s] %s run initiated." % (time.ctime(), prog_name)
else:
    sys.exit("\nUsage: python %s <in.vcf> <in.dbsnp.bad>\n" % prog_name)
# fi


def retrieve_bitarray_data(bad_file):
    bad = dict()
    byte_dict = marshal.load( open(bad_file, "rb") )
    for chrom in byte_dict:
        bad[chrom] = bitarray.bitarray()
        bad[chrom].frombytes( byte_dict[chrom] )
    # for chrom end
    return bad


############
##  MAIN  ##
############

bad_db = retrieve_bitarray_data(in_bad_db)

for line in open(in_vcf, "r"):
    if line.startswith('#'):
        print line.strip()
        continue
    field = line.strip().split('\t')
    chrom = field[0]
    one_pos = int(field[1])
    zero_pos = one_pos - 1

    ## Pedigree filter
    if (bad_db[chrom][zero_pos]):
        continue
    
    print line.strip()
# for line end
