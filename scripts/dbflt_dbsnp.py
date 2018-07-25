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
    print("[%s] %s run initiated." % (time.ctime(), prog_name), file=sys.stderr)
else:
    sys.exit("\nUsage: python %s <in.vcf> <in.db> <maf.cut>\n" % prog_name)
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
        db_alts = result[4]
        if db_pos.count('_') or db_alts.count('<') or db_alts.count('>'):
            iter_cnt -= 1
            continue

        # Inclusion citeria
        db_info = result[7]
        flag_CDA = db_info.count("CDA") # mut in clin diagnostic assay
        flag_PM = db_info.count("PM;") # pubmed mut
        flag_MUT = db_info.count("MUT") # low freq cited mut

        # Exclusion criteria
        flag_G5 = db_info.count("G5") # common SNP
        if flag_G5: # filter out
            print_flag = False
            continue

        db_af_src = re.search("CAF=(.+?);", db_info)
        if db_af_src:
            mafs = [float(x) for x in db_af_src.groups()[0].split(',')
                    if x != '.']
            if len(mafs) > 0:
                max_maf = max(mafs[1:]) if len(mafs)>1 else mafs[0]
            # Exclude by MAF 
            if max_maf > 0.05: # above 0.05 is benign anyways...
                print_flag = False
                continue
            elif max_maf < maf_cut: # DB MAF < MAF cutff
                continue
            else: # exclude if var MAF > MAF cut
                if flag_CDA or flag_PM or flag_MUT: # don't filter out
                    continue
                else:
                    print_flag = False
                    continue
            # fi
        else:
            continue
        # fi
    # for result end
    if print_flag:
        print(line.rstrip())
# for line end
