#!/usr/bin/python

############
##  INIT  ##
############

import os, sys

if len(sys.argv) == 3:
    info = sys.argv[1]
    data_dir = sys.argv[2] 
else:
    prog_name = sys.argv[0].split('/')[-1]
    sys.exit("\nUsage: python %s <in.info.txt> <data.dir>\n" % prog_name)
# fi


#################
##  FUNCTIONS  ## --> or modules, whatever you prefer :)
#################

# Make running options
def get_flt_opt( inh_type ):
    if inh_type == "AD": # autosomal dominant. e.g. Aa x AA --> Aa
        inh_flt = "inherited"
        pat_flt = "het1"
        fat_flt = "het1" if fat_aff=="y" else "ref"
        mot_flt = "het1" if mot_aff=="y" else "ref"
    elif inh_type == "XD": # X-linked dominant. e.g. XY x X'X --> X'X
        inh_flt = "inherited"
        pat_flt = "het1_hom"
        fat_flt = "hom" if fat_aff=="y" else "ref"
        mot_flt = "het1" if mot_aff=="y" else "ref"
    elif inh_type == "AR": # autosomal recessive. e.g. Aa x Aa --> aa
        inh_flt = "inherited"
        pat_flt = "hom"
        fat_flt = "hom" if fat_aff=="y" else "het1"
        mot_flt = "hom" if mot_aff=="y" else "het1"
    elif inh_type == "XR": # X-linked recessive. e.g. X'Y x X'X --> X'X'
        inh_flt = "inherited"
        pat_flt = "hom"
        fat_flt = "hom"
        assert fat_aff=="y", "ERROR: Father should be affected for %s" % inh_type
        mot_flt = "hom" if mot_aff=="y" else "het1"
    elif inh_type == "DN": # de novo. e.g. AA x AA --> Aa
        inh_flt = "de_novo"
        pat_flt = "not_ref"
        fat_flt = "ref"
        mot_flt = "ref"
        assert (fat_aff=="n" and mot_aff=="n"), "ERROR: Parents should NOT be affected for %s" % inh_type
    elif inh_type == "AC": # autosomal compound heterozygous. e.g. A1A2,a1A2 x A1A2,A1a2 --> a1A2,A1a2
        inh_flt = "inherited"
        pat_flt = "het1"
        fat_flt = "all"
        mot_flt = "all"
    elif inh_type == "XC": # X-linked heterozygous. e.g. X1x2,Y x X1X2,x1X2 --> X1x2,x1X2
        inh_flt = "inherited"
        pat_flt = "het1"
        fat_flt = "ref_hom"
        assert fat_aff=="y", "ERROR: Father should be affected for %s" % inh_type
        mot_flt = "all"
    elif inh_type == "YL": # Y-linked. XY' x XX --> XY'
        inh_flt = "inherited"
        pat_flt = "hom"
        fat_flt = "hom"
        assert fat_aff=="y", "ERROR: Father should be affected for %s" % inh_type
        mot_flt = "ref"
    elif inh_type == "MT": # mitochondrial. MMMm x MMMmmm --> mm
        inh_flt = "inherited"
        pat_flt = "not_ref"
        fat_flt = "all"
        mot_flt = "all"
    else:
        sys.exit("ERROR: Unsupported inheritance type (%s)" % inh_type)
    # fi

    return inh_flt, pat_flt, fat_flt, mot_flt
# fed


# Set GMAF cutoff
def get_maf_cut( stringency ):
    if stringency == "4":
        maf_cut = "0"
    elif stringency == "3":
        maf_cut = "0005"
    elif stringency == "2":
        maf_cut = "001"
    elif stringency == "1":
        maf_cut = "010"
    else:
        sys.exit("ERROR: Unsupported stringency (%s)" % stringency)
    # fi
    return maf_cut
# fed


############
##  MAIN  ##
############

# Get path for data
abs_path = os.path.abspath("%s" % sys.argv[0])
abs_dir = abs_path.rsplit('/',1)[0]

# Parse info file
info_file = open(info, "r")
line = True # init
while line != "":
    line = info_file.readline()
    if line.startswith("# Case ID"):
        next_line = info_file.readline()
        while next_line.startswith('#'):
            next_line = info_file.readline()
            continue
        psym = next_line.split('#')[0].strip()
        continue
    elif line.startswith("# Sample ID"):
        next_line = info_file.readline()
        while next_line.startswith('#'):
            next_line = info_file.readline()
            continue
        pat_id, fat_id, mot_id = next_line.split('#')[0].strip().split('/')
        continue
    elif line.startswith("# Affected?"):
        next_line = info_file.readline()
        while next_line.startswith('#'):
            next_line = info_file.readline()
            continue
        pat_aff, fat_aff, mot_aff = next_line.split('#')[0].strip().split('/')
        continue
    elif line.startswith("# Inheritance type"):
        next_line = info_file.readline()
        while next_line.startswith('#'):
            next_line = info_file.readline()
            continue
        inh_types = next_line.split('#')[0].strip().split(',')
        continue
    elif line.startswith("# Filter stringency"):
        next_line = info_file.readline()
        while next_line.startswith('#'):
            next_line = info_file.readline()
            continue
        stringencies = next_line.split('#')[0].strip().split(',')
        continue
    else:
        continue
    # fi
# while end
info_file.close()

# Run code for each INHERITANCE_TYPE and STRINGENCY
for inh_type in inh_types:
    inh_flt, pat_flt, fat_flt, mot_flt = get_flt_opt( inh_type )

    # Check sample IDs
    assert pat_id != "n", "ERROR: Sample ID should be present."
    fat_id = "-" if fat_id=="n" else fat_id
    mot_id = "-" if mot_id=="n" else mot_id
    
    fat_flt = "-" if fat_id=="-" else fat_flt
    mot_flt = "-" if mot_id=="-" else mot_flt
    
    ch_flag = "1" if inh_type.count("compound heterozygote") else "0"
    
    for stringency in stringencies:
        maf_cut = get_maf_cut(stringency)
        # Make command
        cmd = ("/bin/bash %s/scripts/run_var_flt.sh %s %s %s %s %s %s %s %s %s %s %s %s"
               % (abs_dir, psym, inh_flt, ch_flag, maf_cut, pat_id, pat_flt, fat_id, fat_flt, mot_id, mot_flt, inh_type, data_dir))
        #print cmd
        os.system(cmd)
    # for str end
# for inh end
