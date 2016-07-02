# Created: December 2nd 2014
# Last update: July 1st 2016
# Author: Seongmin Choi <seongmin.choi@raregenomics.org>

[ $# -ne 4 ] && { echo -e "\nUsage: $0 <data.dir> <case.symbol> <out.bam.symbol> <id.file>\n"; exit 1; }

abs_path=`readlink -e $0`
scripts_dir=`dirname $abs_path`
rainbow_dir=${scripts_dir%%\/scripts}

data_dir=$1 # data directory, where fastq dir is present
[ ! -d $data_dir ] && { echo "ERROR: $data_dir does not exist."; exit 1; }
data_dir=`readlink -e ${data_dir}`
psym=$2 # case symbol
asym=$3 # output bam symbol (e.g. rg --> sample.rg.bam)
id_file=$4 # includes all the sample IDs for the case
[ ! -f $id_file ] && { echo "ERROR: $id_file does not exist."; exit 1; }

rc=$rainbow_dir/rc # rc file
[ ! -f $rc ] && { echo "ERROR: $rc does not exist." 1>&2; exit 1; }
source $rc # import rc

rglb=$psym
rgpl=illumina
rgpu=$psym

picard=$rainbow_dir/bin/$picard_pf/AddOrReplaceReadGroups.jar # picard run file
[ ! -f $picard ] && { echo "ERROR: $picard does not exist."; exit 1; }

in_dir=$data_dir/bwa_output/$psym # bam file dir
[ ! -d $in_dir ] && { echo "ERROR: $in_dir does not exist."; exit 1; }

out_dir=$data_dir/rg_sort_bam/$psym # output bam dir
[ ! -d $out_dir ] && { mkdir $out_dir; }

log_dir=$out_dir/log # log dir
[ ! -d $log_dir ] && { mkdir $log_dir; }

for id in `cat $id_file`; do
    in_bam=${id}.bam
    [ ! -f $in_dir/$in_bam ] && { echo "ERROR: $in_dir/$in_bam does not exist."; exit 1; }
    out_bam=${id}.${asym}.bam
    log_file=${id}.${asym}.log
    cmd="java -Xmx4g -Djava.io.tmpdir=${out_dir}/tmp -jar $picard"
    cmd="$cmd I=$in_dir/$in_bam"
    cmd="$cmd O=$out_dir/$out_bam"
    cmd="$cmd RGLB=${rglb}"
    cmd="$cmd RGPL=${rgpl}"
    cmd="$cmd RGPU=${rgpu}"
    cmd="$cmd RGSM=${id}"
    cmd="$cmd VALIDATION_STRINGENCY=LENIENT"
    cmd="$cmd SORT_ORDER=coordinate"
    cmd="$cmd CREATE_INDEX=true"
    cmd="$cmd > $log_dir/$log_file 2>&1"

    eval "$cmd"
done

