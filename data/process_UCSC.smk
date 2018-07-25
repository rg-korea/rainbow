import os, re

configfile: config['configfile']
db_name = "UCSC"
wk_dir = os.path.split(os.path.abspath(config['configfile']))[0]
db_dir = os.path.join(wk_dir, config['database_dir'])
db_subdir = os.path.join(db_dir, db_name)

refFlat = os.path.join(db_subdir, 'refFlat.txt.gz')
simpleRepeat = os.path.join(db_subdir, 'simpleRepeat.txt.gz')
refFlat_bed = re.sub(r'\.txt\.gz$', '.bed.gz', refFlat)
refFlat_tbi = re.sub(r'\.txt\.gz$', '.bed.gz.tbi', refFlat)
simpleRepeat_bed = re.sub(r'\.txt\.gz$', '.bed.gz', simpleRepeat)
simpleRepeat_tbi = re.sub(r'\.txt\.gz$', '.bed.gz.tbi', simpleRepeat)

flag_file = os.path.join(db_subdir, 'COMPLETE')

rule all:
    input:
        refFlat_bed,
        refFlat_tbi,
        simpleRepeat_bed,
        simpleRepeat_tbi
    output:
        flag_file
    shell:
        'touch {output}'

rule file_check:
    input:
        refFlat,
        simpleRepeat
    shell:
        'echo "All files present: {input}"'

rule proc_refFlat:
    input:
        refFlat
    output:
        refFlat_bed
    shell:
        'zcat {input} | cut -f3,5,6 | sort -k2,3 -n | bedtools sort | uniq | bgzip -c > {output}'

rule tabix_refFlat:
    input:
        refFlat_bed
    output:
        refFlat_tbi
    shell:
        'tabix -p vcf {input}'

rule proc_simpleRepeat:
    input:
        simpleRepeat
    output:
        simpleRepeat_bed
    shell:
        'zcat {input} | cut -f2,3,4 | bgzip -c > {output}'

rule tabix_simpleRepeat:
    input:
        simpleRepeat_bed
    output:
        simpleRepeat_tbi
    shell:
        'tabix -p vcf {input}'


