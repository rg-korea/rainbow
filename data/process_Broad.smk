import os, re

configfile: config['configfile']
db_name = "Broad"
wk_dir = os.path.split(os.path.abspath(config['configfile']))[0]
db_dir = os.path.join(wk_dir, config['database_dir'])
db_subdir = os.path.join(db_dir, db_name)
db_files = expand(os.path.join(db_subdir, '{fname}'), fname=config['files_'+db_name])

ref_files = [x.rsplit('.gz', 1)[0] for x in db_files 
    if x.endswith('.gz')]
fasta = [x for x in ref_files if (x.endswith('.fasta') or x.endswith('.fa'))][0]
bwa_indices = [fasta+'.amb', fasta+'.ann', fasta+'.bwt', fasta+'.pac', fasta+'.sa']

flag_file = os.path.join(db_subdir, 'COMPLETE')

rule all:
    input:
        bwa_indices
    output:
        flag_file
    shell: 
        "touch {output}"

rule file_check:
    input:
        db_files
    shell:
        'echo "All files present: {input}"'

rule unzip_files:
    input:
        db_files
    output:
        ref_files
    shell:
        'gzip -d {input}'

rule bwa_index:
    input:
        fasta
    output:
        bwa_indices
    shell:
        'bwa index -a bwtsw {input}'

