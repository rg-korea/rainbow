[ $# -ne 2 ] && { echo -e "\nUsage: $0 <project symbol> <in.vep>\n"; exit 1; }
psym=$1
in_vep=$2

echo "[`date`] $0 run initiated." 1>&2

cmd="cat $in_vep"
cmd="$cmd | grep NM_"
cmd="$cmd | egrep"
cmd="$cmd \"transcript_ablation"
cmd="$cmd|splice_acceptor_variant"
cmd="$cmd|splice_donor_variant"
cmd="$cmd|stop_gained"
cmd="$cmd|frameshift_variant"
cmd="$cmd|stop_lost"
cmd="$cmd|start_lost"
cmd="$cmd|transcript_amplification"
cmd="$cmd|inframe_insertion"
cmd="$cmd|inframe_deletion"
cmd="$cmd|missense_variant"
cmd="$cmd|protein_altering_variant"
cmd="$cmd|incomplete_terminal_codon_variant"
cmd="$cmd|coding_sequence_variant\""

[ -f $in_vep ] && { eval $cmd; }
[ ! -f $in_vep ] && { echo ""; }
