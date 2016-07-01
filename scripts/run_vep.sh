# Last update: July 1st 2016
# Author: Seongmin Choi <seongmin.choi@raregenomics.org>

perl=/usr/bin/perl # do not change version from v5.18.2
[ ! -f $perl ] && { echo "ERROR: $perl does not exist."; exit 1; }
bash=/bin/bash
[ ! -f $bash ] && { echo "ERROR: $bash does not exist."; exit 1; }
vep_dir=/BiO/rgi/bin/ensembl-tools-release-78/scripts/variant_effect_predictor
[ ! -d $vep_dir ] && { echo "ERROR: $vep_dir does not exist."; exit 1; }
vep=variant_effect_predictor.pl
[ ! -f $vep_dir/$vep ] && { echo "ERROR: $vep_dir/$vep does not exist."; exit 1; }

cache_dir=/BiO/rgi/data/db/vep
[ ! -d $cache_dir ] && { echo "ERROR: $cache_dir does not exist."; exit 1; }

[ $# -ne 2 ] && { echo -e "\nUsage: $0 <project symbol> <patient ID>\n"; exit 1; }
psym=$1
pat_id=$2

in_dir=/BiO/rgi/data/var_flt/$psym
[ ! -d $in_dir ] && { echo "ERROR: $in_dir does not exist."; exit 1; }

out_dir=$in_dir
[ ! -d $out_dir ] && { echo "LOG: $out_dir does not exist. Creating one..."; mkdir $out_dir; }

log_dir=$out_dir/log
[ ! -d $log_dir ] && { echo "LOG: $log_dir does not exist. Creating one..."; mkdir $log_dir; }

in_vcf=$in_dir/${pat_id}.fin.vcf
[ ! -f $in_vcf ] && { echo "ERROR: $in_vcf does not exist."; exit 1; }

out_vep=$out_dir/${pat_id}.vep
log_file=$log_dir/${pat_id}.vep.log

cmd="$perl $vep_dir/$vep"
cmd="$cmd --cache"
cmd="$cmd --offline"
cmd="$cmd --sift b"
cmd="$cmd --polyphen b"
cmd="$cmd --symbol"
cmd="$cmd --refseq"
cmd="$cmd --force_overwrite"
cmd="$cmd --dir $cache_dir"
cmd="$cmd -i $in_vcf"
cmd="$cmd -o $out_vep"
cmd="$cmd > $log_file 2>&1"

echo $cmd
