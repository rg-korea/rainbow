# Input VCF, input BAD (bitarray dictionary) to select exonic variants.
# BAD format example: bad["chr12"] = bitarray("0000000001111111111000000 ...")
# where the 1's are exon regions made from a Captured_Exon_Regions.bed file.

############
##  INIT  ##
############

import sys, bitarray, marshal, time

prog_name = sys.argv[0].split('/')[-1]
if len(sys.argv) == 3:
    in_vcf = sys.argv[1]
    in_bad = sys.argv[2]
    print >> sys.stderr, "[%s] %s run initiated." % (time.ctime(), prog_name)
else:
    sys.exit("\nUsage: python %s <in.vcf> <in.bad>\n" % prog_name)
# fi


############
##  MAIN  ##
############

def retrieve_bitarray_data(bad_file): # from marshal dictionary to bitarray dictionary
    bad = dict()
    byte_dict = marshal.load( open(bad_file, "rb") )
    for chrom in byte_dict:
        bad[chrom] = bitarray.bitarray()
        bad[chrom].frombytes( byte_dict[chrom] )
    # for chrom end
    return bad


bad_db = retrieve_bitarray_data(in_bad)

for line in open(in_vcf, "r"):
    if line.startswith('#'):
        print line.strip()
        continue

    field = line.strip().split('\t')
    chrom = field[0]
    one_pos = int(field[1])
    zero_pos = one_pos - 1

    ## Mask filter; if chr:pos is in exon pass, else continue
    if (bad_db[chrom][zero_pos]):
        pass
    else:
        continue
    
    print line.strip()
# for line end
