import os

configfile: config['configfile']
db_name = "dbSNP"
wk_dir = os.path.split(os.path.abspath(config['configfile']))[0]
db_dir = os.path.join(wk_dir, config['database_dir'])
db_subdir = os.path.join(db_dir, db_name)
db_files = expand(os.path.join(db_subdir, '{fname}'), fname=config['files_'+db_name])

flag_file = os.path.join(db_subdir, 'COMPLETE')

rule all:
    input:
        db_files
    output:
        flag_file
    shell:
        'touch {output}'

rule file_check:
    input:
        db_files
    shell:
        'echo "All files present: {input}"'
