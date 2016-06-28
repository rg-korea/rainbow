#!/bin/bash

# Settings
perl=/usr/bin/perl
[ ! -f $perl ] && { echo "ERROR: $perl does not exist." 1>&2; exit 1; }

wkdir=`pwd` # working directory
bin=$wkdir/bin
[ ! -d $bin ] && { mkdir $bin; }

rc=./rc # rc file
[ ! -f $rc ] && { echo "ERROR: $rc does not exist." 1>&2; exit 1; }
source $rc # import rc

# Download fasta files from the Broad Institute bundle
[ ! -d $wkdir/data ] && { mkdir $wkdir/data; }
[ ! -d $wkdir/data/index ] && { mkdir $wkdir/data/index; }
for gz_file in ucsc.hg19.dict.gz ucsc.hg19.fasta.fai.gz ucsc.hg19.fasta.gz; do
    wget $broad_ftp/$gz_file -P $wkdir/data/index
    gzip -d $wkdir/data/index/$gz_file
done

# Install python2.7 (user)
wget $python_url -P $bin
[ ! -f $bin/$python_pf.tgz ] && { echo "ERROR: $bin/$python_pf.tgz does not exist." 1>&2; exit 1; }
tar -zxf $bin/$python_pf.tgz -C $bin
cd $bin/$python_pf
[ ! -d $bin/user_python_dir ] && { mkdir $bin/user_python_dir; }
./configure --enable-unicode=ucs4 --prefix=$bin/user_python_dir
make
make install
cd $wkdir
rm $bin/$python_pf.tgz

# Install bitarray (python module)
wget $bitarray_url -P $bin
[ ! -f $bin/$bitarray_pf.tar.gz ] && { echo "ERROR: $bin/$bitarray_pf.tar.gz does not exist." 1>&2; exit 1; }
tar -zxf $bin/$bitarray_pf.tar.gz -C $bin
cd $bin/$bitarray_pf
[ ! -f $bin/$python_pf/python ] && { echo "ERROR: $bin/$python_pf/python does not exist." 1>&2; exit 1; }
$bin/$python_pf/python ./setup.py install
cd $wkdir
rm $bin/$bitarray_pf.tar.gz

# Install bwa (must have "mem" option)
wget $bwa_url -P $bin
[ ! -f $bin/$bwa_pf.tar.bz2 ] && { echo "ERROR: $bin/$bwa_pf.tar.bz2 does not exist." 1>&2; exit 1; }
bzip2 -d $bin/$bwa_pf.tar.bz2
tar -xf $bin/$bwa_pf.tar -C $bin
cd $bin/$bwa_pf
make
[ ! -f $wkdir/data/index/ucsc.hg19.fasta ] && { echo "ERROR: $wkdir/data/index/ucsc.hg19.fasta does not exist." 1>&2; exit 1; }
$bin/$bwa_pf/bwa index -a bwtsw $wkdir/data/index/ucsc.hg19.fasta
cd $wkdir
rm $bin/$bwa_pf.tar

# Install picard-tools-1.96 (must have separate .jar files)
wget $picard_url -P $bin
[ ! -f $bin/$picard_pf.zip ] && { echo "ERROR: $bin/$picard_pf.zip does not exist." 1>&2; exit 1; }
unzip $bin/$picard_pf.zip -d $bin
cd $wkdir
rm $bin/$picard_pf.zip

# Install samtools
wget $samtools_url -P $bin
[ ! -f $bin/$samtools_pf.tar.bz2 ] && { echo "ERROR: $bin/samtools-1.3.tar.bz2 does not exist." 1>&2; exit 1; }
bzip2 -d $bin/$samtools_pf.tar.bz2
tar -xf $bin/$samtools_pf.tar -C $bin
cd $bin/$samtools_pf
make
cd $wkdir
rm $bin/$samtools_pf.tar

# Install Platypus
## TO-DO

# Install htslib
wget $htslib_url -P $bin
bzip2 -d $bin/$htslib_pf.tar.bz2
tar -xf $bin/$htslib_pf.tar -C $bin
cd $bin/$htslib_pf
mkdir $HOME/rainbow
./configure --prefix=$HOME/rainbow
export PATH=$PATH:$HOME/rainbow/bin
make
make install
rm $bin/$htslib_pf.tar

# Install Ensembl Varient Effect Predictor 78
wget $vep_url -P $bin
[ ! -f $bin/$vep_pf.zip ] && { echo "ERROR: $bin/$vep_pf.zip does not exist." 1>&2; exit 1; }
unzip $bin/$vep_pf.zip -d $bin
cd $bin/ensembl-tools-release-$vep_pf/scripts/variant_effect_predictor
[ ! -d $wkdir/data/db ] && { mkdir $wkdir/data/db; }
[ ! -d $wkdir/data/db/vep ] && { mkdir $wkdir/data/db/vep; }
printf 'y\n43' | $perl $bin/ensembl-tools-release-$vep_pf/scripts/variant_effect_predictor/INSTALL.pl -c $wkdir/data/db/vep --SPECIES homo_sapiens_refseq --ASSEMBLY GRCh37
cd $wkdir
rm $bin/$vep_pf.zip

### Download bitarray dictionaries: SNP filter database (in-house)
##[ ! -d $wkdir/data/db/ba_dict ] && { mkdir $wkdir/data/db/ba_dict; }
##wget -r $bitarray_dict_url -P $wkdir/data/db/ba_dict

