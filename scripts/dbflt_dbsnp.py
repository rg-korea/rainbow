# Created: September 12th 2016
# Last update: September 19th 2016
# Author: Seongmin Choi <seongmin.choi@raregenomics.org>

############
##  INIT  ##
############

import sys, tabix, time, re

prog_name = sys.argv[0].split('/')[-1]
if len(sys.argv) == 4:
    in_vcf = sys.argv[1]
    in_db = sys.argv[2]
    maf_cut = float(sys.argv[3])
    print >> sys.stderr, "[%s] %s run initiated." % (time.ctime(), prog_name)
else:
    sys.exit("\nUsage: python %s <in.vcf> <in.db.url> <maf.cut>\n" % prog_name)
# fi

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
    query_db = "%s:%s-%s" % (chrom_id, one_pos, one_pos)
    iter_cnt = 0
    
    try:
        results = db.querys(query_db) # send query
    except tabix.TabixError:
        print line.strip()
        continue

    # Process tabix results
    for result in results:
        iter_cnt += 1
        db_chrom_id = result[0]
        db_pos = result[1]
        db_alts = result[4]
        if db_pos.count('_') or db_alts.count('<') or db_alts.count('>'):
            iter_cnt -= 1
            continue

        db_info = result[7]
        db_af_src = re.search("CAF=([\d\.,-]+);", db_info)

        # Inclusion citeria
        flag_CDA = db_info.count("CDA") # mut in clin diagnostic assay
        flag_PM = db_info.count("PM;") # precious mut
        flag_MUT = db_info.count("MUT") # low freq cited mut
        if flag_CDA or flag_PM or flag_MUT: # don't filter out
            print line.strip()
            continue

        # Exclusion criteria
        flag_G5 = db_info.count("G5") # common SNP
        if flag_G5: # filter out
            continue

        if db_af_src:
            mafs = [float(x) for x in db_af_src.groups()[0].split(',')
                    if x != '.']
            max_maf = max(mafs)
            # Exclude by MAF 
            if max_maf < maf_cut: # DB MAF < MAF cutff
                print line.strip()
                continue
            else: # exclude if var MAF > MAF cut
                continue
            # fi
        else:
            print line.strip()
            continue
    # for result end

    # No result found in DB
    if iter_cnt == 0:
        print line.strip()
        continue
    else:
        continue
# for line end
