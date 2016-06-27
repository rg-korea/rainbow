[ $# -ne 1 ] && { echo -e "\nUsage: $0 <in.vcf>\n"; exit 1; }

vcf=$1
[ ! -f $vcf ] && { echo "ERROR: $vcf does not exist."; exit 1; }

grep -v "^#" $vcf | awk '{ printf("%s\t%s\t%s\n", $1, $2-1, $2) }'
