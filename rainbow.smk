import os

configfile: config['configfile']

db_dir = config['database_dir']
case = config["case"]
data_dir = config["data_dir"]
maf_cutoff = str(config["maf_cutoff"])
patient_id = config["patient_id"]
patient_gt = config["patient_gt"]
father_id = config["father_id"]
father_gt = config["father_gt"]
mother_id = config["mother_id"]
mother_gt = config["mother_gt"]
smp_ids = [patient_id, father_id, mother_id]
pedigree = config["pedigree"]
inh_filter = config["inh_filter"]
filter_out_non_omim = config["filter_out_non_omim"]
omim_genemap = config["omim_genemap2_path"] # optional 

ref_fasta = os.path.join(db_dir, "Broad", "ucsc.hg19.fasta")

def get_fastq_path(data_dir, smp_id):
    fq1 = os.path.join(data_dir, smp_id+"_1.fastq.gz")
    fq2 = os.path.join(data_dir, smp_id+"_2.fastq.gz")
    if not os.path.isfile(fq1):
        fq1 = os.path.join(data_dir, smp_id+"_1.fastq")
        if not os.path.isfile(fq1):
            fq1 = ''
    if not os.path.isfile(fq2):
        fq2 = os.path.join(data_dir, smp_id+"_1.fastq")
        if not os.path.isfile(fq2):
            fq2 = ''
    return fq1, fq2

patient_fq1, patient_fq2 = get_fastq_path(data_dir, patient_id)
patient_bam = os.path.join(data_dir, patient_id+".bam")
patient_log = os.path.join(data_dir, "log", patient_id+".bam.log")
patient_addrg_bam = os.path.join(data_dir, patient_id+".addrg.bam")
patient_addrg_log = os.path.join(data_dir, "log", patient_id+".addrg.bam.log")
patient_samtools_idx = os.path.join(data_dir, patient_id+".addrg.bam.bai")
patient_samtools_log = os.path.join(data_dir, "log", patient_id+".addrg.bam.bai.log")

father_fq1, father_fq2 = get_fastq_path(data_dir, father_id)
father_bam = os.path.join(data_dir, father_id+".bam")
father_log = os.path.join(data_dir, "log", father_id+".bam.log")
father_addrg_bam = os.path.join(data_dir, father_id+".addrg.bam")
father_addrg_log = os.path.join(data_dir, "log", father_id+".addrg.bam.log")
father_samtools_idx = os.path.join(data_dir, father_id+".addrg.bam.bai")
father_samtools_log = os.path.join(data_dir, "log", father_id+".addrg.bam.bai.log")

mother_fq1, mother_fq2 = get_fastq_path(data_dir, mother_id)
mother_bam = os.path.join(data_dir, mother_id+".bam")
mother_log = os.path.join(data_dir, "log", mother_id+".bam.log")
mother_addrg_bam = os.path.join(data_dir, mother_id+".addrg.bam")
mother_addrg_log = os.path.join(data_dir, "log", mother_id+".addrg.bam.log")
mother_samtools_idx = os.path.join(data_dir, mother_id+".addrg.bam.bai")
mother_samtools_log = os.path.join(data_dir, "log", mother_id+".addrg.bam.bai.log")

caller_inputs = [patient_addrg_bam]
if os.path.isfile(mother_fq1):
    print(mother_fq1, "is a file")
    caller_inputs.append(mother_addrg_bam)
if os.path.isfile(father_fq1):
    print(father_fq1, "is a file")
    caller_inputs.append(father_addrg_bam)

case_vcf = os.path.join(data_dir, case+".vcf")
case_log = os.path.join(data_dir, "log", case+".log")
snpeff_vcf = os.path.join(data_dir, case+".snpEff.vcf")
rflt_vcf = os.path.join(data_dir, case+".rflt.vcf")
dbflt_dbSNP = case+".dbflt_dbSNP.MAF"+maf_cutoff+".vcf"
dbflt_dbSNP = os.path.join(data_dir, dbflt_dbSNP)
dbflt_1000G = case+".dbflt_1000G.MAF"+maf_cutoff+".vcf"
dbflt_1000G = os.path.join(data_dir, dbflt_1000G)
dbflt_GoNL = case+".dbflt_GoNL.MAF"+maf_cutoff+".vcf"
dbflt_GoNL = os.path.join(data_dir, dbflt_GoNL)
dbflt_EVS = case+".dbflt_EVS.MAF"+maf_cutoff+".vcf"
dbflt_EVS = os.path.join(data_dir, dbflt_EVS)
dbflt_ExAC = case+".dbflt_ExAC.MAF"+maf_cutoff+".vcf"
dbflt_ExAC = os.path.join(data_dir, dbflt_ExAC)
dbflt_STR = case+".dbflt_STR.MAF"+maf_cutoff+".vcf"
dbflt_STR = os.path.join(data_dir, dbflt_STR)
spec_vcf = os.path.join(data_dir, case+".spec.MAF"+maf_cutoff+".vcf")
pflt_vcf = os.path.join(data_dir, patient_id+".pflt.MAF"+maf_cutoff+"."+pedigree+".vcf")
raw_csv = os.path.join(data_dir, patient_id+".MAF"+maf_cutoff+"."+pedigree+".raw.csv")
omim_annot_csv = os.path.join(data_dir, patient_id+".MAF"+maf_cutoff+"."+pedigree+".omim_annot.csv")

db_dbSNP = os.path.join(db_dir, "dbSNP", config['files_dbSNP'][0])
db_1000G = os.path.join(db_dir, "1000G")
db_GoNL = os.path.join(db_dir, "GoNL")
db_EVS = os.path.join(db_dir, "EVS")
db_ExAC = os.path.join(db_dir, "ExAC", config['files_ExAC'][0])

refFlat = os.path.join(db_dir, "UCSC", config['files_UCSC'][0])
simpleRepeat = os.path.join(db_dir, "UCSC", config['files_UCSC'][1])

rule all:
    input: 
        raw_csv if omim_genemap == '-' else omim_annot_csv

rule bwa_patient:
    input: ref_fasta, patient_fq1, patient_fq2
    output: patient_bam
    log: patient_log
    shell:
        "(bwa mem {input} | samtools view -Sb - > {output}) 2> {log}"

rule addrg_patient:
    input: patient_bam
    output: patient_addrg_bam
    params:
        rgpl="-PL illumina",
        rglb=expand("-LB {patient_id}", patient_id=patient_id),
        rgpu=expand("-PU {patient_id}", patient_id=patient_id),
        rgsm=expand("-SM {patient_id}", patient_id=patient_id),
        tmp_dir=expand("{data_dir}/tmp", data_dir=data_dir)
    log: patient_addrg_log
    shell:
        "gatk --java-options '-Xmx4g -Djava.io.tmpdir={params.tmp_dir}' "
        "AddOrReplaceReadGroups -I {input} -O {output} --CREATE_INDEX=true "
        "{params.rgpl} {params.rglb} {params.rgpu} {params.rgsm} "
        "-SO coordinate > {log} 2>&1"

rule samtools_index_patient:
    input: patient_addrg_bam
    output: patient_samtools_idx
    log: patient_samtools_log
    shell:
        "samtools index -b {input} {output} > {log} 2>&1"

if father_id != '-':
    rule bwa_father:
        input: ref_fasta, father_fq1, father_fq2
        output: father_bam
        log: father_log
        shell:
            "(bwa mem {input} | samtools view -Sb - > {output}) 2> {log}"

    rule addrg_father:
        input: father_bam
        output: father_addrg_bam
        params:
            rgpl="-PL illumina",
            rglb=expand("-LB {father_id}", father_id=father_id),
            rgpu=expand("-PU {father_id}", father_id=father_id),
            rgsm=expand("-SM {father_id}", father_id=father_id),
            tmp_dir=expand("{data_dir}/tmp", data_dir=data_dir)
        log: father_addrg_log
        shell:
            "gatk --java-options '-Xmx4g -Djava.io.tmpdir={params.tmp_dir}' "
            "AddOrReplaceReadGroups -I {input} -O {output} --CREATE_INDEX=true "
            "{params.rgpl} {params.rglb} {params.rgpu} {params.rgsm} "
            "-SO coordinate > {log} 2>&1"

    rule samtools_index_father:
        input: father_addrg_bam
        output: father_samtools_idx
        log: father_samtools_log
        shell:
            "samtools index -b {input} {output} > {log} 2>&1"

if mother_id != '-':
    rule bwa_mother:
        input: ref_fasta, mother_fq1, mother_fq2
        output: mother_bam
        log: mother_log
        shell:
            "(bwa mem {input} | samtools view -Sb - > {output}) 2> {log}"

    rule addrg_mother:
        input: mother_bam
        output: mother_addrg_bam
        params:
            rgpl="-PL illumina",
            rglb=expand("-LB {mother_id}", mother_id=mother_id),
            rgpu=expand("-PU {mother_id}", mother_id=mother_id),
            rgsm=expand("-SM {mother_id}", mother_id=mother_id),
            tmp_dir=expand("{data_dir}/tmp", data_dir=data_dir)
        log: mother_addrg_log
        shell:
            "gatk --java-options '-Xmx4g -Djava.io.tmpdir={params.tmp_dir}' "
            "AddOrReplaceReadGroups -I {input} -O {output} --CREATE_INDEX=true "
            "{params.rgpl} {params.rglb} {params.rgpu} {params.rgsm} "
            "-SO coordinate > {log} 2>&1"

    rule samtools_index_mother:
        input: mother_addrg_bam
        output: mother_samtools_idx
        log: mother_samtools_log
        shell:
            "samtools index -b {input} {output} > {log} 2>&1"

rule haplotypecaller:
    input: caller_inputs
    params: 
        ref=ref_fasta,
        inputs=expand("-I {bam}", bam=caller_inputs)
    output: case_vcf
    log: case_log
    shell:
        "gatk --java-options '-Xmx4g -Djava.io.tmpdir=./tmp' HaplotypeCaller {params.inputs} -O {output} -R {params.ref} > {log} 2>&1"
rule snpeff:
    input: case_vcf
    output: snpeff_vcf
    log: case_log
    shell:
        "snpEff -noStats hg19 {input} > {output} 2>> {log}"

rule region_filter:
    input: snpeff_vcf, refFlat
    output: rflt_vcf
    log: case_log
    shell:
        "python scripts/dbflt_refflat.py "
        "{input} > {output} 2>> {log}"

rule dbSNP_filter:
    input: rflt_vcf, db_dbSNP
    output: dbflt_dbSNP
    params: maf_cutoff
    log: case_log
    shell:
        "python scripts/dbflt_dbsnp.py "
        "{input} {params} > {output} 2>> {log}"

##@##rule oneKG_filter:
##@##    input: dbflt_dbSNP, db_1000G
##@##    output: dbflt_1000G
##@##    params: maf_cutoff
##@##    log: case_log
##@##    shell:
##@##        "python scripts/dbflt_1000g.py "
##@##        "{input} {params} > {output} 2>> {log}"
##@##
##@##rule GoNL_filter:
##@##    input: dbflt_1000G, db_GoNL
##@##    output: dbflt_GoNL
##@##    params: maf_cutoff
##@##    log: case_log
##@##    shell:
##@##        "python scripts/dbflt_gonl.py "
##@##        "{input} {params} > {output} 2>> {log}"
##@##
##@##rule EVS_filter:
##@##    input: dbflt_GoNL, db_EVS
##@##    output: dbflt_EVS
##@##    params: maf_cutoff
##@##    log: case_log
##@##    shell:
##@##        "python scripts/dbflt_evs.py "
##@##        "{input} {params} > {output} 2>> {log}"
##@##
##@##rule ExAC_filter:
##@##    input: dbflt_EVS, db_ExAC
##@##    output: dbflt_ExAC
##@##    params: maf_cutoff
##@##    log: case_log
##@##    shell:
##@##        "python scripts/dbflt_exac.py "
##@##        "{input} {params} > {output} 2>> {log}"
##@##
##@##rule STR_filter:
##@##    input: dbflt_ExAC, simpleRepeat
##@##    output: dbflt_STR
##@##    log: case_log
##@##    shell:
##@##        "python scripts/dbflt_str.py "
##@##        "{input} > {output} 2>> {log}"

rule specify_variants:
    #input: dbflt_STR
    input: dbflt_dbSNP
    output: spec_vcf
    params: father_id, mother_id
    log: case_log
    shell:
        "python scripts/specify_variants.py "
        "{input} {params} "
        "> {output} 2>> {log}"

rule pedigree_filter:
    input: spec_vcf
    output: pflt_vcf
    params:
        str(inh_filter), pedigree,
        patient_id + ':' + patient_gt,
        father_id + ':' + father_gt,
        mother_id + ':' + mother_gt
    log: case_log + '.pflt.log'
    shell:
        "python scripts/pedigree_flt.py "
        "{input} {params} "
        "> {output} 2>> {log}"

rule make_table:
    input: pflt_vcf
    output: raw_csv
    params: patient_id
    log: case_log
    shell:
        "python scripts/make_table.py "
        "{input} {params} "
        "> {output} 2>> {log}"

if omim_genemap != '-' and os.path.exists(omim_genemap):
    rule add_omim_annotation:
        input: raw_csv, omim_genemap
        output: omim_annot_csv
        params: filter_out_non_omim
        log: case_log
        shell:
            "python scripts/add_omim_to_csv.py "
            "{input} {params} "
            "> {output} 2>> {log}"

