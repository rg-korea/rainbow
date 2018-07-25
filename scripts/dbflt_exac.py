############
##  INIT  ##
############

from __future__ import print_function
import sys, tabix, time, re

prog_name = sys.argv[0].split('/')[-1]
if len(sys.argv) == 4:
    in_vcf = sys.argv[1]
    in_db = sys.argv[2]
    maf_cut = float(sys.argv[3])
    print(sys.stderr, "[%s] %s run initiated." % (time.ctime(), prog_name), file=sys.stderr)
else:
    sys.exit("\nUsage: python %s <in.vcf> <in.db.url> <maf.cut>\n" % prog_name)
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
    query_db = "%s:%s-%s" % (chrom_id, one_pos, one_pos)
    iter_cnt = 0
    
    try:
        results = db.querys(query_db) # send query
    except tabix.TabixError:
        print(line.strip())
        continue

    # Process tabix results
    print_flag = True
    for result in results:
        iter_cnt += 1
        db_chrom_id = result[0]
        db_pos = result[1]
        db_res = result[3]
        db_alts = result[4].split(',')
        db_chr_pos = "%s:%s" % (db_chrom_id, db_pos)

        db_info = result[7]
        db_af_src = re.search(";AF=(.+?);", db_info)
        if db_af_src:
            mafs = [float(x) for x in db_af_src.groups()[0].split(',')]
            max_maf = max(mafs)
            # Filter by MAF 
            if max_maf > 0.05:
                print_flag = False
                continue
            if max_maf < maf_cut: # DB MAF < MAF cutff
                continue
            else:
                print_flag = False
                continue
            # fi
        else:
            print_flag = False
            continue
    # for result end

    # No result found in DB
    if print_flag:
        print(line.strip())
# for line end
