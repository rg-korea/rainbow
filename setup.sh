#!/bin/bash
# Created: December 23nd 2015
# Last update: July 8th 2016
# Author: Seongmin Choi <seongmin.choi@raregenomics.org>

# Settings
perl=/usr/bin/perl
[ ! -f $perl ] && { echo "ERROR: $perl does not exist." 1>&2; exit 1; }

wkdir=`pwd` # working directory
bin=$wkdir/bin
[ ! -d $bin ] && { mkdir $bin; }

pm=$HOME/perl5 # perl modules config
[ ! -d $pm ] && { mkdir $pm; }
export PERL5LIB=$pm/share/perl5:$PERL5LIB
export PERL5LIB=$pm/lib64/perl5:$PERL5LIB
export PERL5LIB=$pm/lib/perl5:$PERL5LIB

rc=./rc # rc file
[ ! -f $rc ] && { echo "ERROR: $rc does not exist." 1>&2; exit 1; }
source $rc # import rc

setup_log=$wkdir/setup.log
echo "[`date`] Initiating setup.sh." > $setup_log

# Functions
function write_error {
    log_msg="[`date`] ERROR: $1"
    echo $log_msg >> $2
    echo $log_msg 1>&2
    exit 1
}

function write_log {
    log_msg="[`date`] LOG: $1"
    echo $log_msg >> $2
    echo $log_msg 1>&2
    return 0
}

function check_dir {
    if [ ! -d $1 ]; then {
        write_error "$1 does not exist." $setup_log
        exit 1
    }; fi
    return 0
}

function check_file {
    if [ ! -f $1 ]; then {
        write_error "$1 does not exist." $setup_log
        exit 1
    }; fi
    return 0
}

# Download fasta files from the Broad Institute bundle
[ ! -d $wkdir/data ] && { mkdir $wkdir/data; }
[ ! -d $wkdir/data/index ] && { mkdir $wkdir/data/index; }
for gz_file in ucsc.hg19.dict.gz ucsc.hg19.fasta.fai.gz ucsc.hg19.fasta.gz; do
    idx_file=${gz_file%%\.gz}
    if [ -f $wkdir/data/index/$idx_file ]; then 
        write_log "$idx_file already exists." $setup_log
    elif [ -f $wkdir/data/index/$gz_file ]; then {
        write_log "$gz_file already exists." $setup_log &&
        gzip - d $wkdir/data/index/$gz_file &&
        write_log "Completed index file unzip." $setup_log ||
        write_error "Index file not unzipped." $setup_log
    } 
    else {
        wget $broad_ftp/$gz_file -P $wkdir/data/index &&
        write_log "Completed fasta download." $setup_log ||
        write_error "Index file not downloaded." $setup_log
        gzip -d $wkdir/data/index/$gz_file &&
        write_log "Completed index file unzip." $setup_log ||
        write_error "Index file not unzipped." $setup_log
    }; fi
done

# Install python2.7 (user)
res=`$bin/$python_pf/python --version 2>&1 | grep "Python 2.7" | wc -l`
if [ $res -ne 0 ]; then
    write_log "Python already exists." $setup_log
else {
    wget $python_url -P $bin; check_file $bin/$python_pf.tgz &&
    tar -zxf $bin/$python_pf.tgz -C $bin &&
    cd $bin/$python_pf &&
    [ ! -d $bin/user_python_dir ] && { mkdir $bin/user_python_dir; } &&
    ./configure --enable-unicode=ucs4 --prefix=$bin/user_python_dir &&
    make &&
    make install &&
    cd $wkdir &&
    rm $bin/$python_pf.tgz &&
    write_log "Completed python install." $setup_log || write_error "Python not installed." $setup_log
}; fi

# Install bitarray (python module)
res=`$bin/$python_pf/python -c "import bitarray; bitarray.test()" 2>&1 | grep OK | wc -l`
if [ $res -ne 0 ]; then
    write_log "bitarray already exists." $setup_log
else {
    wget $bitarray_url -P $bin &&
    check_file $bin/$bitarray_pf.tar.gz &&
    tar -zxf $bin/$bitarray_pf.tar.gz -C $bin &&
    cd $bin/$bitarray_pf &&
    check_file $bin/$python_pf/python &&
    $bin/$python_pf/python ./setup.py install &&
    cd $wkdir &&
    rm $bin/$bitarray_pf.tar.gz &&
    write_log "Completed bitarray install." $setup_log || write_error "bitarray not installed." $setup_log
}; fi

# Install bwa (must have "mem" option)
res=`$bin/$bwa_pf/bwa 2>&1 | grep BWA-MEM | wc -l`
if [ $res -ne 0 ]; then
    write_log "BWA already exists." $setup_log
else {
    wget $bwa_url -P $bin &&
    check_file $bin/$bwa_pf.tar.bz2 &&
    bzip2 -d $bin/$bwa_pf.tar.bz2 &&
    tar -xf $bin/$bwa_pf.tar -C $bin &&
    cd $bin/$bwa_pf &&
    make &&
    check_file $wkdir/data/index/ucsc.hg19.fasta &&
    bin/$bwa_pf/bwa index -a bwtsw $wkdir/data/index/ucsc.hg19.fasta &&
    cd $wkdir &&
    rm $bin/$bwa_pf.tar &&
    write_log "Completed BWA install." $setup_log || write_error "BWA not installed." $setup_log
}; fi

# Install picard-tools-1.96 (must have separate .jar files)
res=`java -jar $bin/$picard_pf/AddOrReplaceReadGroups.jar --version 2>&1 | wc -l`
if [ $res -ne 0 ]; then
    write_log "Picard already exists." $setup_log
else {
    wget $picard_url -P $bin &&
    check_file $bin/$picard_pf.zip &&
    unzip $bin/$picard_pf.zip -d $bin &&
    cd $wkdir &&
    rm $bin/$picard_pf.zip &&
    write_log "Completed Picard install." $setup_log || write_error "Picard not installed." $setup_log
}; fi

# Install samtools
res=`$bin/$samtools_pf/samtools --version 2>&1 | grep samtools | wc -l`
if [ $res -ne 0 ]; then
    write_log "SAMtools already exists." $setup_log
else {
    wget $samtools_url -P $bin &&
    check_file $bin/$samtools_pf.tar.bz2 &&
    bzip2 -d $bin/$samtools_pf.tar.bz2 &&
    tar -xf $bin/$samtools_pf.tar -C $bin &&
    cd $bin/$samtools_pf &&
    make &&
    cd $wkdir &&
    rm $bin/$samtools_pf.tar &&
    write_log "Completed SAMtools install."  $setup_log || write_error "SAMtools not installed."  $setup_log
}; fi

# Install Platypus
res=`$bin/$python_pf/python $bin/$platypus_pf/Platypus.py callVariants --help | grep Usage | wc -l`
if [ $res -ne 0 ]; then
    write_log "Platypus already exists." $setup_log
else {
    check_file $bin/$platypus_pf.tar.gz &&
    tar -zxf $bin/$platypus_pf.tar.gz -C $bin &&
    cd $bin/$platypus_pf &&
    ./buildPlatypus.sh &&
    cd $wkdir &&
    write_log "Completed Platypus install." $setup_log || write_error "Platypus not installed." $setup_log
}; fi

# Install Recursive.pm
res=`ls $bin/$fileref_pf | wc -l` 
if [ $res -ne 0 ]; then
    write_log "Recursive.pm already exists." $setup_log
else {
    wget $filerec_url -P $bin &&
    check_file $bin/$filerec_pf.tar.gz &&
    tar -zxf $bin/$filerec_pf.tar.gz -C $bin &&
    cd $bin/$filerec_pf &&
    perl Makefile.PL PREFIX=$pm &&
    make &&
    make install &&
    cd $wkdir &&
    write_log "Completed Recursive.pm install." $setup_log || write_error "Recursive.pm not installed." $setup_log
}; fi

# Install htslib
res=`ls $bin/$htslib_pf | wc -l` 
if [ $res -ne 0 ]; then
    write_log "HTSlib already exists." $setup_log
else {
    wget $htslib_url -P $bin &&
    check_file $bin/$htslib_pf.tar.bz2 &&
    bzip2 -d $bin/$htslib_pf.tar.bz2 &&
    tar -xf $bin/$htslib_pf.tar -C $bin &&
    cd $bin/$htslib_pf &&
    mkdir $HOME/rainbow &&
    check_dir $HOME/rainbow &&
    ./configure --prefix=$HOME/rainbow &&
    export PATH=$PATH:$HOME/rainbow/bin &&
    make &&
    make install &&
    rm $bin/$htslib_pf.tar &&
    cd $wkdir &&
    write_log "Completed HTSlib install." $setup_log || write_error "HTSlib not installed." $setup_log
}; fi

# Install DBI
res=`ls $bin/$dbi_pf | wc -l` 
if [ $res -ne 0 ]; then
    write_log "DBI.pm already exists." $setup_log
else {
    wget $dbi_url -P $bin &&
    check_file $bin/$dbi_pf.tar.gz &&
    tar -xzf $bin/$dbi_pf.tar.gz -C $bin &&
    cd $bin/$dbi_pf &&
    perl Makefile.PL INSTALL_BASE=$pm &&
    make &&
    make install &&
    cd $wkdir &&
    write_log "Completed DBI.pm install." $setup_log || write_error "DBI.pm not installed." $setup_log
}; fi

# Install BioPerl
res=`ls $bin/$bioperl_pf | wc -l`
if [ $res -ne 0 ]; then
    write_log "BioPerl already exists." $setup_log
else {
    wget $bioperl_url -P $bin &&
    check_file $bin/$bioperl_pf.tar.gz &&
    tar -zxf $bin/$bioperl_pf.tar.gz -C $bin &&
    cd $bin/$bioperl_pf &&
    printf 'anyyyyyyyyyyyyyyyyyyyyyyyyy' | perl Build.PL PREFIX=$pm LIB=$pm &&
    yes y | ./Build installdeps &&
    ./Build install &&
    cd $wkdir &&
    write_log "Completed BioPerl install." $setup_log || write_error "BioPerl not installed." $setup_log
}; fi

# Install Ensembl Varient Effect Predictor
[ ! -d $wkdir/data/db ] && { mkdir $wkdir/data/db; } &&
check_dir $wkdir/data/db &&
[ ! -d $wkdir/data/db/vep ] && { mkdir $wkdir/data/db/vep; } &&
check_dir $wkdir/data/db/vep &&
refseq_db_dir=$wkdir/data/db/vep/homo_sapiens_refseq
res=`ls $refseq_db_dir | grep GRCh37 | wc -l`
if [ $res -ne 0 ]; then
    write_log "VEP already exists." $setup_log
else {
    ensembl_dir=$bin/ensembl-tools-release-$vep_pf &&
    vep_dir=$ensembl_dir/scripts/variant_effect_predictor &&
    [ -d $ensembl_dir ] && { rm -rf $ensembl_dir; } &&
    [ -d $refseq_db_dir ] && { rm -rf $refseq_db_dir; } &&
    wget $vep_url -P $bin &&
    check_file $bin/$vep_pf.zip &&
    unzip $bin/$vep_pf.zip -d $bin &&
    export PERL5LIB=$PERL5LIB:$vep_dir &&
    cd $vep_dir &&
    $perl $bin/ensembl-tools-release-$vep_pf/scripts/variant_effect_predictor/INSTALL.pl -a ac -c $wkdir/data/db/vep -s homo_sapiens_refseq -y GRCh37 &&
    cd $wkdir &&
    rm $bin/$vep_pf.zip &&
    write_log "Completed VEP install." $setup_log || write_error "VEP not installed." $setup_log
}; fi

### Download bitarray dictionaries: SNP filter database (in-house)
##[ ! -d $wkdir/data/db/ba_dict ] && { mkdir $wkdir/data/db/ba_dict; }
##wget -r $bitarray_dict_url -P $wkdir/data/db/ba_dict

