############
##  INIT  ##
############

from __future__ import print_function
import sys, tabix, time, re

prog_name = sys.argv[0].split('/')[-1]
if len(sys.argv) == 3:
    in_vcf = sys.argv[1]
    in_db = sys.argv[2] # bed.gz
    print("[%s] %s run initiated." % (time.ctime(), prog_name), file=sys.stderr)
else:
    sys.exit("\nUsage: python %s <in.vcf> <in.db.bed.gz>\n" % prog_name)
# fi

# Init tabix
db = tabix.open(in_db)

# Proc VCF
for line in open(in_vcf, "r"):
    if line.startswith('#'):
        print(line.strip())
        continue
    field = line.strip().split('\t')
    chrom = field[0]
    chrom_id = chrom.replace("chr", '')
    chrom_id = 'M' if chrom_id == "MT" else chrom_id
    one_pos = int(field[1])
    chr_pos = "%s:%s" % (chrom_id, one_pos)
    ref = field[3]
    alts = field[4].split(',')
    query_db = "chr%s:%s-%s" % (chrom_id, one_pos, one_pos)
    try:
        results = db.querys(query_db) # send query
        iter_cnt = sum(1 for _ in results)
    except tabix.TabixError:
        print(line.strip())
        continue

    # If at least 1 gene region present
    if iter_cnt > 0:
        print(line.strip())
        continue
    elif iter_cnt == 0:
        continue
    else:
        sys.exit("ERROR: iter_cnt = %s" % iter_cnt)
# for line end
