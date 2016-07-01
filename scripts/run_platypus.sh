# Last update: July 1st 2016
# Author: Seongmin Choi <seongmin.choi@raregenomics.org>

[ $# -ne 4 ] && { echo -e "\nUsage: $0 <data.dir> <case.symbol> <in.bam.symbol> <id.file>\n"; exit 1; }

abs_path=`readlink -e $0`
scripts_dir=`dirname $abs_path`
rainbow_dir=${scripts_dir%%\/scripts}

data_dir=$1 # data directory, where fastq dir is present
[ ! -d $data_dir ] && { echo "ERROR: $data_dir does not exist."; exit 1; }
data_dir=`readlink -e ${data_dir}`
psym=$2 # case symbol; lca, cas
asym=$3
id_file=$4 # includes all the sample IDs for the case
[ ! -f $id_file ] && { echo "ERROR: $id_file does not exist."; exit 1; }

rc=$rainbow_dir/rc # rc file
[ ! -f $rc ] && { echo "ERROR: $rc does not exist." 1>&2; exit 1; }
source $rc # import rc

python=$rainbow_dir/bin/$python_pf/python
[ ! -f $python ] && { echo "ERROR: $python does not exist."; exit 1; }
program=$rainbow_dir/bin/$platypus_pf/Platypus.py
[ ! -f $program ] && { echo "ERROR: $program does not exist."; exit 1; }

ref=$rainbow_dir/data/index/ucsc.hg19.fasta # hg19 fasta
[ ! -f $ref ] && { echo "ERROR: $ref does not exist."; exit 1; }

in_dir=$data_dir/rg_sort_bam/$psym # sample fastq file dir
[ ! -d $in_dir ] && { echo "ERROR: $in_dir does not exist."; exit 1; }

out_dir=$data_dir/platypus/$psym # output dir
[ ! -d $out_dir ] && { mkdir $out_dir; }

log_dir=$out_dir/log # log dir
[ ! -d $log_dir ] && { mkdir $log_dir; }

thread=1 # number of threads for multi-processing
min_flank=0

bamfiles=""
for id in `cat $id_file`; do
    in_bam=${id}.${asym}.bam
    [ ! -f $in_dir/$in_bam ] && { echo "ERROR: $in_dir/$in_bam does not exist."; exit 1; }
    bamfiles="${bamfiles},${in_dir}/${in_bam}"
done
bamfiles=${bamfiles##,}

out_file=${psym}.vcf
log_file=${out_file}.log

cmd="$python"
cmd="$cmd $program"
cmd="$cmd callVariants"
cmd="$cmd --refFile $ref"
cmd="$cmd --bamFiles $bamfiles"
cmd="$cmd --output $out_dir/$out_file"
cmd="$cmd --nCPU $thread"
cmd="$cmd --minFlank $min_flank"
cmd="$cmd --logFileName $log_dir/$log_file"
cmd="$cmd > $log_dir/$log_file 2>&1"

eval "$cmd"
