# Created: September 12th 2016
# Last update: September 19th 2016
# Author: Seongmin Choi <seongmin.choi@raregenomics.org>

############
##  INIT  ##
############

import sys, tabix, time, re

prog_name = sys.argv[0].split('/')[-1]
if len(sys.argv) == 3:
    in_vcf = sys.argv[1]
    in_db = sys.argv[2] # bed.gz
    print >> sys.stderr, "[%s] %s run initiated." % (time.ctime(), prog_name)
else:
    sys.exit("\nUsage: python %s <in.vcf> <in.db.bed.gz>\n" % prog_name)
# fi


## Module
def has_indel(ref, alts):
    flag_has_indel = False
    ref_len = len(ref)

    for alt in alts:
        if ref_len == len(alt) == 1:
            continue
        else:
            flag_has_indel = True
        # fi
    # for end 
    return flag_has_indel
# fed

# Init tabix
db = tabix.open(in_db)

# Proc VCF
for line in open(in_vcf, "r"):
    if line.startswith('#'):
        print line.strip()
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
    flag_has_indel = has_indel(ref, alts)
    if not flag_has_indel: # if var has no indel
        print line.strip()
        continue

    try:
        results = db.querys(query_db) # send query
        iter_cnt = sum(1 for _ in results)
    except tabix.TabixError:
        print line.strip()
        continue

    # If at least 1 STR present
    if iter_cnt > 0:
        continue
    elif iter_cnt == 0:
        print line.strip()
        continue
    else:
        sys.exit("ERROR: iter_cnt = %s" % iter_cnt)
# for line end
