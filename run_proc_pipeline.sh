#!/bin/bash
# Last update: July 1st 2016
# Author: Seongmin Choi <seongmin.choi@raregenomics.org>

[ $# -ne 4 ] && { echo -e "\nUsage: $0 <data.dir> <case.symbol> <end.form> <id.file>\n" 1>&2; exit 1; }

# Settings
abs_path=`readlink -e $0`
rainbow_dir=`dirname $abs_path`

bash=/bin/bash
[ ! -f $bash ] && { echo "ERROR: $bash does not exist." 1>&2; exit 1; }

input_dir=$1 # data directory
[ ! -d $data_dir ] && { echo "ERROR: $input_dir does not exist." 1>&2; exit 1; }
data_dir=`readlink -e ${input_dir}`

psym=$2 # case symbol

endform=$3

id_file=$4 # includes all the sample IDs for the case

rc=$rainbow_dir/rc # rc file
[ ! -f $rc ] && { echo "ERROR: $rc does not exist." 1>&2; exit 1; }
source $rc # import rc

bwa=$rainbow_dir/bin/$bwa_pf/bwa # bwa run file (must have "mem" option)

samtools=$rainbow_dir/bin/$samtools_pf/samtools # samtools run file
[ ! -f $samtools ] && { echo "ERROR: $samtools does not exist." 1>&2; exit 1; }

picard=$rainbow_dir/bin/$picard_pf/AddOrReplaceReadGroups.jar # picard run file
[ ! -f $picard ] && { echo "ERROR: $picard does not exist." 1>&2; exit 1; }

bwa_dir=$data_dir/bwa_output
[ ! -d $bwa_dir ] && { mkdir $bwa_dir; }

rg_dir=$data_dir/rg_sort_bam
[ ! -d $rg_dir ] && { mkdir $rg_dir; }

vcf_dir=$data_dir/platypus
[ ! -d $vcf_dir ] && { mkdir $vcf_dir; }

vcf_flt_dir=$data_dir/var_flt
[ ! -d $vcf_flt_dir ] && { mkdir $vcf_flt_dir; }

result_dir=$data_dir/result
[ ! -d $result_dir ] && { mkdir $result_dir; }


# Make run commands
cmd1="$bash $rainbow_dir/scripts/run_bwa_mem.sh $data_dir $psym $endform $id_file"

cmd2="$bash $rainbow_dir/scripts/run_addrg.sh $data_dir $psym rg $id_file"

cmd3="$bash $rainbow_dir/scripts/run_ln_bai_bambai.sh $data_dir/rg_sort_bam/$psym; sleep 5"

cmd4="$bash $rainbow_dir/scripts/run_platypus.sh $data_dir $psym rg $id_file"

cmds="$cmd1; $cmd2; $cmd3; $cmd4"
eval $cmds # run commands linearly
