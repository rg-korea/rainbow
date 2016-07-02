# Created: October 30th 2014
# Last update: July 1st 2016
# Author: Seongmin Choi <seongmin.choi@raregenomics.org>

[ $# -ne 4 ] && { echo -e "\nUsage: $0 <working.dir> <in.dir> <out.dir>\n"; exit 1; }
rc=$wkdir/rc # rc file
[ ! -f $rc ] && { echo "ERROR: $rc does not exist." 1>&2; exit 1; }
source $rc # import rc

samtools=$wkdir/$samtools_pf/samtools # samtools run file

in_dir=$1 # bam file dir
out_dir=$2 # output index dir

for bam in `ls $in_dir | egrep "\.bam$"`; do
    cmd="$samtools index"
    cmd="$cmd $in_dir/$bam"
    cmd="$cmd $out_dir/${bam%%bam}bai"
    eval "$cmd"
done
