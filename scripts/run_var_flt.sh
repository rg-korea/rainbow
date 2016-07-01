# Last update: July 1st 2016
# Author: Seongmin Choi <seongmin.choi@raregenomics.org>

# Settings
perl=/usr/bin/perl # do not change version from v5.18.2
[ ! -f $perl ] && { echo "ERROR: $perl does not exist."; exit 1; }

bash=/bin/bash
[ ! -f $bash ] && { echo "ERROR: $bash does not exist."; exit 1; }

abs_path=`readlink -e $0`
scripts_dir=`dirname $abs_path`
[ ! -d $scripts_dir ] && { echo "ERROR: $scripts_dir does not exist." 1>&2; exit 1; }
rainbow_dir=${scripts_dir%%\/scripts}
[ ! -d $rainbow_dir ] && { echo "ERROR: $rainbow_dir does not exist." 1>&2; exit 1; }

bin=$rainbow_dir/bin
[ ! -d $bin ] && { mkdir $bin; }

rc=./rc # rc file
[ ! -f $rc ] && { echo "ERROR: $rc does not exist." 1>&2; exit 1; }
source $rc # import rc


# Input
[ $# -ne 12 ] && { echo -e "\nUsage: $0 <case.symbol> <inh.flt> <comp.het?> <AF.cut> <pat.id> <pat.flt> <fat.id> <fat.flt> <mot.id> <mot.flt> <inh.type> <data.dir>\n"; exit 1; }
psym=$1 # case symbol; lca, cas
inh_flt=$2
ch_flag=$3
af_cut=$4
pat_id=$5
pat_flt=$6
fat_id=$7
fat_flt=$8
mot_id=$9
mot_flt=${10}
inh_type=${11}
data_dir=${12}


# Check dependencies
python=$rainbow_dir/bin/$python_pf/python
[ ! -f $python ] && { echo "ERROR: $python does not exist."; exit 1; }
vep_dir=$rainbow_dir/bin/ensembl-tools-release-$vep_pf/scripts/variant_effect_predictor
[ ! -d $vep_dir ] && { echo "ERROR: $vep_dir does not exist."; exit 1; }
vep=variant_effect_predictor.pl
[ ! -f $vep_dir/$vep ] && { echo "ERROR: $vep_dir/$vep does not exist."; exit 1; }
vep_opts="--cache --offline --sift b --polyphen b --symbol --refseq --force_overwrite"

var_flt_dir=$data_dir/var_flt
raw_vcf_dir=$data_dir/platypus/$psym
flt_vcf_dir=$data_dir/var_flt/$psym
result_dir=$data_dir/result
ba_dict_dir=$rainbow_dir/data/db/ba_dict
cache_dir=$rainbow_dir/data/db/vep
[ ! -d $var_flt_dir ] && { echo "LOG: $var_flt_dir does not exist. Creating directory..."; mkdir $var_flt_dir; }
[ ! -d $raw_vcf_dir ] && { echo "ERROR: $raw_vcf_dir does not exist."; exit 1; }
[ ! -d $flt_vcf_dir ] && { echo "LOG: $flt_vcf_dir does not exist. Creating directory..."; mkdir $flt_vcf_dir; }
[ ! -d $flt_vcf_dir ] && { echo "LOG: $result_dir does not exist. Creating directory..."; mkdir $flt_vcf_dir; }
[ ! -d $ba_dict_dir ] && { echo "ERROR: $ba_dict_dir does not exist."; exit 1; }
[ ! -d $cache_dir ] && { echo "ERROR: $cache_dir does not exist."; exit 1; }

db1=$ba_dict_dir/dbsnp_138.hg19.caf${af_cut}.vcf.bad
db2=$ba_dict_dir/1000G.hg19.maf${af_cut}.vcf.bad
db3=$ba_dict_dir/gonl.r5.hg19.maf${af_cut}.vcf.bad
db4=$ba_dict_dir/snp142Common.flt.bed.bad
db5=$ba_dict_dir/EVS.hg19.maf${af_cut}.vcf.bad
db6=$ba_dict_dir/ExAC.r0.3.sites.vep.hg19.maf${af_cut}.vcf.bad
[ ! -f $db1 ] && { echo "ERROR: $db1 does not exist."; exit 1; }
[ ! -f $db2 ] && { echo "ERROR: $db2 does not exist."; exit 1; }
[ ! -f $db3 ] && { echo "ERROR: $db3 does not exist."; exit 1; }
[ ! -f $db4 ] && { echo "ERROR: $db4 does not exist."; exit 1; }
[ ! -f $db5 ] && { echo "ERROR: $db5 does not exist."; exit 1; }
[ ! -f $db6 ] && { echo "ERROR: $db6 does not exist."; exit 1; }

region_bad=$ba_dict_dir/hg19_refGene_ext7.bad # ref_gene region
[ ! -f $region_bad ] && { echo "ERROR: $region_bad does not exist."; exit 1; }

raw_vcf=$raw_vcf_dir/$psym.vcf
[ ! -f $raw_vcf ] && { echo "ERROR: $raw_vcf does not exist."; exit 1; }


# Set output files
reg_flt_vcf=$flt_vcf_dir/${psym}.rflt.vcf # after filtering variants not in exome-captured regions
spec_var_vcf=$flt_vcf_dir/${psym}.spec.vcf # after specifying variants
var_flt_vcf=$flt_vcf_dir/${pat_id}.${inh_type}_${af_cut}.vflt.vcf # after filtering low quality variants
ped_flt_vcf=$flt_vcf_dir/${pat_id}.${inh_type}_${af_cut}.pflt.vcf # after filtering variants not segregated according to the pedigree
db1_flt_vcf=$flt_vcf_dir/${pat_id}.${inh_type}_${af_cut}.d1flt.vcf # after filtering common variants with DB1
db2_flt_vcf=$flt_vcf_dir/${pat_id}.${inh_type}_${af_cut}.d2flt.vcf # after filtering common variants with DB2
db3_flt_vcf=$flt_vcf_dir/${pat_id}.${inh_type}_${af_cut}.d3flt.vcf # after filtering common variants with DB3
db4_flt_vcf=$flt_vcf_dir/${pat_id}.${inh_type}_${af_cut}.d4flt.vcf # after filtering common variants with DB4
db5_flt_vcf=$flt_vcf_dir/${pat_id}.${inh_type}_${af_cut}.d5flt.vcf # after filtering common variants with DB5
db6_flt_vcf=$flt_vcf_dir/${pat_id}.${inh_type}_${af_cut}.d6flt.vcf # after filtering common variants with DB6
cmp_het_vcf=$flt_vcf_dir/${pat_id}.${inh_type}_${af_cut}.cpflt.vcf # after selecting compound hetero variants
fin_vcf=$flt_vcf_dir/${pat_id}.${inh_type}_${af_cut}.fin.vcf # variants of the patient
raw_vep=$flt_vcf_dir/${pat_id}.${inh_type}_${af_cut}.vep
flt_vep=$flt_vcf_dir/${pat_id}.${inh_type}_${af_cut}.flt.vep
cmp_het_vep=$flt_vcf_dir/${pat_id}.${inh_type}_${af_cut}.ch.vep 
cmp_het_flt_vep=$flt_vcf_dir/${pat_id}.${inh_type}_${af_cut}.cp.flt.vep
fin_output=$result_dir/${pat_id}.${inh_type}_${af_cut}.result.txt


# Set commands
cmd1="$python $scripts_dir/get_target_var.py $raw_vcf $region_bad > $reg_flt_vcf" # filter variants not in the exome-captured regions
cmd2="$python $scripts_dir/specify_variants.py $reg_flt_vcf $fat_id $mot_id > $spec_var_vcf" # filter variants according to the given pedigree
cmd3="$python $scripts_dir/platypus_call_flt.py $spec_var_vcf $pat_id > $var_flt_vcf" # filter variants according to calling quality
cmd4="$python $scripts_dir/pedigree_flt.py $var_flt_vcf $inh_flt $inh_type $pat_id:$pat_flt $fat_id:$fat_flt $mot_id:$mot_flt > $ped_flt_vcf" # filter variants according to the given pedigree
cmd5="$python $scripts_dir/dbsnp_flt.py $ped_flt_vcf $db1 > $db1_flt_vcf" # filter variants according to common SNP data
cmd6="$python $scripts_dir/dbsnp_flt.py $db1_flt_vcf $db2 > $db2_flt_vcf" # filter variants according to common SNP data
cmd7="$python $scripts_dir/dbsnp_flt.py $db2_flt_vcf $db3 > $db3_flt_vcf" # filter variants according to common SNP data
cmd8="$python $scripts_dir/dbsnp_flt.py $db3_flt_vcf $db4 > $db4_flt_vcf" # filter variants according to common SNP data
cmd9="$python $scripts_dir/dbsnp_flt.py $db4_flt_vcf $db5 > $db5_flt_vcf" # filter variants according to common SNP data
cmd10="$python $scripts_dir/dbsnp_flt.py $db5_flt_vcf $db6 > $db6_flt_vcf" # filter variants according to common SNP data
cmd11="$python $scripts_dir/extract_patient_vcf.py $db6_flt_vcf $pat_id > $fin_vcf" # final VCF
cmd12="echo \"[`date`] VEP run initiated.\" 1>&2; $perl $vep_dir/$vep $vep_opts --dir $cache_dir -i $fin_vcf -o $raw_vep > /dev/null 2>&1" # run VEP
cmd13="$bash $scripts_dir/vep_flt.sh $psym $raw_vep > $flt_vep" # filter VEP output
cmd14="$python $scripts_dir/arrange_vep.py $flt_vep $fin_vcf $pat_id > $fin_output" # make final output

cmds="$cmd1; $cmd2; $cmd3; $cmd4; $cmd5; $cmd6; $cmd7; $cmd8; $cmd9; $cmd10; $cmd11; $cmd12; $cmd13; $cmd14" # commands for non-compound-heteros

cmdF="$cmd1; $cmd2; $cmd3; $cmd4; $cmd5; $cmd6; $cmd7; $cmd8; $cmd9; $cmd10; $cmd11; $cmd12; $cmd13" # commands for compound heteros - part 1
cmdC1="$python $scripts_dir/compound_hetero_flt.py $db6_flt_vcf $flt_vep $pat_id $fat_id $mot_id > $cmp_het_vcf" # take compound hetero variant only
cmdC2="$python $scripts_dir/extract_patient_vcf.py $cmp_het_vcf $pat_id > $fin_vcf" # final VCF
cmdC3="echo \"[`date`] VEP run initiated.\" 1>&2; $perl $vep_dir/$vep $vep_opts --dir $cache_dir -i $fin_vcf -o $cmp_het_vep > /dev/null 2>&1" # run VEP agian, for compound heteros
cmdC4="$bash $scripts_dir/vep_flt.sh $psym $cmp_het_vep > $cmp_het_flt_vep" # filter VEP output
cmdC5="$python $scripts_dir/arrange_vep.py $cmp_het_flt_vep $fin_vcf $pat_id > $fin_output" # make final output
[ $ch_flag -ne 0 ] && { cmds="$cmdF; $cmdC1; $cmdC2; $cmdC3; $cmdC4; $cmdC5"; } # commands for compound heteros - part 2

# Run
cmds="$cmds"
eval $cmds
