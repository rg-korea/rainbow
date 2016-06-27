[ $# -ne 1 ] && { echo -e "\nUsage: $0 <in.dir>\n"; exit 1; }

in_dir=$1 # bai file dir
[ ! -d $in_dir ] && { echo "ERROR: $in_dir does not exist."; exit 1; }
in_dir=`readlink -e ${in_dir}`
echo $in_dir

for bai in `ls $in_dir | egrep "\.rg.bai$"`; do
    in_bai=$in_dir/$bai
    out_bai=$in_dir/${bai%bai}bam.bai
    [ ! -f $in_bai ] && { echo "ERROR: $in_bai does not exist." 1>&2; continue; } 
    [ -h $out_bai ] && { echo "LOG: $out_bai already exists." 1>&2; continue; } 

    cmd="ln -s $in_bai"
    cmd="$cmd $out_bai"
    eval $cmd
done

