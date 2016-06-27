thread=1 # number of threads for multi-processing

[ $# -ne 4 ] && { echo -e "\nUsage: $0 <data.dir> <case.symbol> <se/pe> <id.file>\n"; exit 1; }

abs_path=`readlink -e $0`
scripts_dir=`dirname $abs_path`
rainbow_dir=${scripts_dir%%\/scripts}

data_dir=$1 # data directory, where fastq dir is present
[ ! -d $data_dir ] && { echo "ERROR: $data_dir does not exist."; exit 1; }
data_dir=`readlink -e ${data_dir}`
psym=$2 # case symbol; lca, cas
endform=$3 # se, pe
id_file=$4 # includes all the sample IDs for the case.
[ ! -f $id_file ] && { echo "ERROR: $id_file does not exist."; exit 1; }

rc=$rainbow_dir/rc # rc file
[ ! -f $rc ] && { echo "ERROR: $rc does not exist." 1>&2; exit 1; }
source $rc # import rc

bwa=$rainbow_dir/bin/$bwa_pf/bwa # bwa run file
samtools=$rainbow_dir/bin/$samtools_pf/samtools # samtools run file
[ ! -f $bwa ] && { echo "ERROR: $bwa does not exist."; exit 1; }
[ ! -f $samtools ] && { echo "ERROR: $samtools does not exist."; exit 1; }

ref=$rainbow_dir/data/index/ucsc.hg19.fasta # hg19 fasta
[ ! -f $ref ] && { echo "ERROR: $ref does not exist."; exit 1; }

in_dir=$data_dir/fastq/$psym # sample fastq file dir
[ ! -d $in_dir ] && { echo "ERROR: $in_dir does not exist."; exit 1; }

out_dir=$data_dir/bwa_output/$psym # output dir
[ ! -d $out_dir ] && { mkdir $out_dir; }

log_dir=$out_dir/log # log dir
[ ! -d $log_dir ] && { mkdir $log_dir; }

if [ $endform = "se" ]; then # single-end
    for id in `cat $id_file`; do
        in_fastq=${id}.fastq.gz
        [ ! -f $in_dir/$in_fastq ] && { echo "ERROR: $in_fastq does not exist."; exit 1; }
        out_bam=${id}.bam
        log_file=${id}.bwa.log
        cmd="$bwa"
        cmd="$cmd mem"
        cmd="$cmd -t $thread"
        cmd="$cmd $ref"
        cmd="$cmd $in_dir/$in_fastq"
        cmd="$cmd 2> $log_dir/$log_file"
        cmd="$cmd | $samtools view -bS - > $out_dir/$out_bam"
        cmd="$cmd 2>> $log_dir/$log_file"
    
        eval "$cmd"
    done
fi

if [ $endform = "pe" ]; then # paired-end
    for id in `cat $id_file`; do
        in_fastq_1=${id}_1.fastq.gz
        in_fastq_2=${id}_2.fastq.gz
        [ ! -f $in_dir/$in_fastq_1 ] && { echo "ERROR: $in_fastq_1 does not exist."; exit 1; }
        [ ! -f $in_dir/$in_fastq_2 ] && { echo "ERROR: $in_fastq_2 does not exist."; exit 1; }
        out_bam=${id}.bam
        log_file=${id}.bwa.log
        cmd="$bwa"
        cmd="$cmd mem"
        cmd="$cmd -t $thread"
        cmd="$cmd $ref"
        cmd="$cmd $in_dir/$in_fastq_1"
        cmd="$cmd $in_dir/$in_fastq_2"
        cmd="$cmd 2> $log_dir/$log_file"
        cmd="$cmd | $samtools view -bS - > $out_dir/$out_bam"
        cmd="$cmd 2>> $log_dir/$log_file"
    
        eval "$cmd"
    done
fi
