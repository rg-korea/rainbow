############
##  INIT  ##
############

from __future__ import print_function
import sys, tabix, time, re, glob

prog_name = sys.argv[0].split('/')[-1]
if len(sys.argv) == 4:
    in_vcf = sys.argv[1]
    in_db_dir = sys.argv[2]
    maf_cut = float(sys.argv[3])
    print(sys.stderr, "[%s] %s run initiated." % (time.ctime(), prog_name), file=sys.stderr)
else:
    sys.exit("\nUsage: python %s <in.vcf> <in.GoNL.vcf.gz.dir> <maf.cut>\n" % prog_name)
# fi

# Init tabix
dbs = {}
for chrom_id in [str(x) for x in range(1, 23)]: # no X/Y/mito var in GoNL
    file_to_glob = "%s/gonl.chr%s.*.vcf.gz" % (in_db_dir, chrom_id)
    db_file = glob.glob(file_to_glob)[0]
    dbs[chrom_id] = tabix.open(db_file)
#db = tabix.open(in_db)

# Proc VCF
for line in open(in_vcf, "r"):
    flag_printed = False
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

    if chrom_id in dbs:
        db = dbs[chrom_id]
        try:
            results = db.querys(query_db) # send query
        except tabix.TabixError:
            print(line.strip())
            continue
        # yrt
    else: # chr non-match; no filtering possible
        print(line.strip())
        continue
    
    # Process tabix results
    print_flag = True
    for result in results:
        iter_cnt += 1
        db_info = result[7]
        db_af_src = re.search("AF=(.+?);", db_info)

        if db_af_src:
            mafs = [float(x) for x in db_af_src.groups()[0].split(',')
                    if x != '.']
            max_maf = max(mafs)
            # Exclude by MAF 
            if max_maf > 0.05:
                print_flag = False
                continue
            elif max_maf < maf_cut: #and not #flag_printed: # DB MAF < MAF cutff
                continue
            else: # exclude if var MAF > MAF cut
                print_flag = False
                continue
        else:
            continue
    # for result end

    # No result found in DB
    if print_flag:
        print(line.rstrip())
# for line end
