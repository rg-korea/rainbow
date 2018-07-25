import os
import sys
import argparse
import configparser
import yaml

# Change java memory options in snpEff
sed_cmd = ("sed -i --follow-symlinks 's/default_jvm_mem_opts="
        '"-Xms512m -Xmx1g"/default_jvm_mem_opts="-Xmx4g"'
        "/g' $(which snpEff)")

# Make running options
def get_flt_opt(config):
    pedigree = config.pedigree
    fat_aff = config.father_aff
    mot_aff = config.mother_aff
    if pedigree == "AD": # autosomal dominant. e.g. Aa x AA --> Aa
        pat_flt = "het1"
        fat_flt = "het1" if fat_aff=="yes" else "ref"
        mot_flt = "het1" if mot_aff=="yes" else "ref"
    elif pedigree == "XD": # X-linked dominant. e.g. XY x X'X --> X'X
        pat_flt = "het1_hom"
        fat_flt = "hom" if fat_aff=="yes" else "ref"
        mot_flt = "het1" if mot_aff=="yes" else "ref"
    elif pedigree == "AH": # autosomal homozygous. e.g. Aa x Aa --> aa
        pat_flt = "hom"
        fat_flt = "hom" if fat_aff=="yes" else "het1"
        mot_flt = "hom" if mot_aff=="yes" else "het1"
    elif pedigree == "XH": # X-linked homozygous. e.g. XY x X'X --> X'Y
        pat_flt = "hom"
        fat_flt = "hom" if fat_aff=="yes" else "ref"
        mot_flt = "hom" if mot_aff=="yes" else "het1"
    elif pedigree == "DN": # de novo. e.g. AA x AA --> Aa
        pat_flt = "not_ref"
        fat_flt = "ref"
        mot_flt = "ref"
        assert (fat_aff!="yes" and mot_aff!="yes"), "ERROR: Parents should NOT be affected for %s" % pedigree
    elif pedigree == "AC": # autosomal compound heterozygous. e.g. A1A2,a1A2 x A1A2,A1a2 --> a1A2,A1a2
        pat_flt = "het1"
        fat_flt = "all"
        mot_flt = "all"
    elif pedigree == "XC": # X-linked heterozygous. e.g. X1x2,Y x X1X2,x1X2 --> X1x2,x1X2
        pat_flt = "het1"
        fat_flt = "ref_hom"
        assert fat_aff=="yes", "ERROR: Father should be affected for %s" % pedigree
        mot_flt = "all"
    elif pedigree == "YL": # Y-linked. XY' x XX --> XY'
        pat_flt = "hom"
        fat_flt = "hom"
        assert fat_aff=="yes", "ERROR: Father should be affected for %s" % pedigree
        mot_flt = "ref"
    elif pedigree == "MT": # mitochondrial. MMMm x MMMmmm --> mm
        pat_flt = "not_ref"
        fat_flt = "all"
        mot_flt = "all"
    else:
        sys.exit("ERROR: Unsupported inheritance type (%s)" % pedigree)
    # fi
    return pat_flt, fat_flt, mot_flt
# fed

def get_cmd(config):
    params = config['PARAMETERS']
    config.case = params['case'].strip()
    config.data_dir = params['data_dir'].strip()
    assert os.path.isdir(config.data_dir), "ERROR: {} does not exist.".format(data_dir)
    config.patient_id = params['patient_id'].strip()
    
    if 'father_id' in params:
        config.father_id = params['father_id'].strip()
        config.father_aff = params['father_affected'].strip()
    else:
        config.father_id = '-'
        config.father_aff = 'no'
    if config.father_id == 'none':
        config.father_id = '-'
    
    if 'mother_id' in params:
        config.mother_id = params['mother_id'].strip()
        config.mother_aff = params['mother_affected'].strip()
    else:
        config.mother_id = '-'
        config.mother_aff = 'no'
    if config.mother_id == 'none':
        config.mother_id = '-'
    
    config.maf_cutoff = params['maf_cutoff'].strip()
    config.pedigree = params['pedigree'].strip()
    config.pedigrees = ['DN', 'AD', 'AH', 'XD', 'XH', 'AC', 'XC', 'YL', 'MT']
    if config.pedigree == "DN":
        config.inh_filter = '2' # de novo
    elif config.pedigree.endswith('C'):
        config.inh_filter = '1' # compound hetero
    elif config.pedigree in config.pedigrees:
        config.inh_filter = '0' # inherited
    else:
        sys.exit("ERROR: {} not in {}.".format(config.pedigree, config.pedigrees))
    
    optionals = config['OPTIONAL']
    if ('omim_genemap2_path' in optionals 
        and os.path.isfile(optionals['omim_genemap2_path'])):
        config.omim_genemap2_path = optionals['omim_genemap2_path']
    else:
        config.omim_genemap2_path = '-'
    config.filter_out_non_omim = optionals['filter_out_non_omim']
    if config.filter_out_non_omim == "yes":
        config.filter_out_non_omim = '1'
    else:
        config.filter_out_non_omim = '0'
    config.database_yaml = optionals['data_yaml']
    assert os.path.isfile(config.database_yaml), "ERROR: {} does not exist.".format(config.database_yaml)
    
    config.pat_flt, config.fat_flt, config.mot_flt = get_flt_opt(config)
    config.mot_flt = '-' if config.mother_id == '-' else config.mot_flt
    config.fat_flt = '-' if config.father_id == '-' else config.fat_flt
    
    config.snakefile = 'rainbow.smk'
    assert os.path.isfile(config.snakefile), "ERROR: {} does not exist.".format(config.snakefile)
    cmd = "snakemake -j -s " + config.snakefile
    cmd += " --config "
    cmd += " configfile=" + config.database_yaml
    cmd += " case=" + config.case
    cmd += " data_dir=" + config.data_dir
    cmd += " maf_cutoff=" + config.maf_cutoff
    cmd += " patient_id=" + config.patient_id
    cmd += " father_id=" + config.father_id
    cmd += " mother_id=" + config.mother_id
    cmd += " patient_gt=" + config.pat_flt
    cmd += " father_gt=" + config.fat_flt
    cmd += " mother_gt=" + config.mot_flt
    cmd += " pedigree=" + config.pedigree
    cmd += " inh_filter=" + config.inh_filter
    cmd += " filter_out_non_omim=" + config.filter_out_non_omim
    cmd += " omim_genemap2_path=" + config.omim_genemap2_path
    return cmd 

if __name__ == "__main__":
    os.system(sed_cmd)
    parser = argparse.ArgumentParser(description="Run rainbow pipeline.")
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('-i', '--ini', required=True,
        help='input a config (.ini)')
    args = parser.parse_args()
    
    config = configparser.ConfigParser()
    assert os.path.isfile(args.ini), "ERROR: "+args.ini+" does not exist."
    config.read(args.ini)
    cmd = get_cmd(config)
   
    print("Running the following snakemake command:\n%s" % cmd)
    os.system(cmd)

