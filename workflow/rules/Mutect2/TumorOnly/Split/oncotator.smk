# A rule to generate a bed from mutect2 vcf, on tumor only 
rule get_variant_bed_tumor_only:
    input:
        Mutect2_vcf = "Mutect2_T/{tsample}_tumor_only_twicefiltered_T.vcf.gz",
        Mutect2_vcf_index = "Mutect2_T/{tsample}_tumor_only_twicefiltered_T.vcf.gz.tbi"
    output:
        BED = "variant_bed_T/{tsample}_tumor_only_T.bed"
    log:
        "logs/variant_bed_T/{tsample}_tumor_only_T.bed.txt"
    params:
        queue = "mediumq",
        vcf2bed = config["APP_VCF2BED"] 
    threads : 1
    resources:
        mem_mb = 5000
    shell:
        'zcat {input.Mutect2_vcf} | python2 {params.vcf2bed} - > {output.BED} 2> {log}'
        
# Run samtools mpileup, on tumor only 
rule samtools_mpileup_tumor_only:
    input:
        BED = "variant_bed_T/{tsample}_tumor_only_T.bed",
        GENOME_REF_FASTA = config["GENOME_FASTA"],
        BAM = "bam/{tsample}.nodup.recal.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{tsample}.recal.bam",
        BAI = "bam/{tsample}.nodup.recal.bam.bai" if config["REMOVE_DUPLICATES"]=="True" else "bam/{tsample}.recal.bam.bai"
    output:
        PILEUP = "pileup_T/{tsample}_tumor_only_T.pileup.gz"
    log:
        "logs/pileup_T/{tsample}_tumor_only_T.pileup.txt"
    params:
        queue = "mediumq",
        samtools = config["APP_SAMTOOLS"]
    threads : 1
    resources:
        mem_mb = 5000
    shell:
        '{params.samtools} mpileup -a -B -l {input.BED} -f {input.GENOME_REF_FASTA} {input.BAM} | gzip - > {output.PILEUP} 2> {log}'
 
## A rule to split mutect2 results in pieces 
rule split_Mutect2_tumor_only:
    input:
        Mutect2_vcf = "Mutect2_T/{tsample}_tumor_only_twicefiltered_T.vcf.gz",
        vcf_index = "Mutect2_T/{tsample}_tumor_only_twicefiltered_T.vcf.gz.tbi",
        interval = config["MUTECT_INTERVAL_DIR"] + "/{interval}.bed"
    output:
        interval_vcf_bcftools = temp("Mutect2_T_oncotator_tmp/{tsample}_tumor_only_T_ON_{interval}_bcftools.vcf.gz"),
        interval_vcf = temp("Mutect2_T_oncotator_tmp/{tsample}_tumor_only_T_ON_{interval}.vcf.gz")
    params:
        queue = "shortq",
        bcftools = config["APP_BCFTOOLS"],
        reformat_script = config["APP_REFORMAT_MUTECT2"]
    log:
        "logs/Mutect2_T_oncotator_tmp/{tsample}_tumor_only_T_ON_{interval}.vcf.log"
    threads : 1
    resources:
        mem_mb = 2000
    shell:
        '{params.bcftools} view -l 9 -R {input.interval} -o {output.interval_vcf_bcftools} {input.Mutect2_vcf} 2> {log}  &&' 
        ' python {params.reformat_script} {output.interval_vcf_bcftools} {output.interval_vcf} 2>> {log}'  
 
# A rule to annotate mutect2 tumor only results with oncotator 
rule oncotator_tumor_only:
    input:
        ONCOTATOR_DB = config["ONCOTATOR_DB"],
        interval_vcf = "Mutect2_T_oncotator_tmp/{tsample}_tumor_only_T_ON_{interval}.vcf.gz"
    output:
        MAF = temp("oncotator_T_tmp/{tsample}_tumor_only_ON_{interval}_annotated_T.TCGAMAF")
    params:
        queue = "mediumq",
        #activate_oncotator = config["ONCOTATOR_ENV"],
        oncotator = config["APP_ONCOTATOR"]
    log:
        "logs/oncotator_T_tmp/{tsample}_ON_{interval}_tumor_only_annotated_T.TCGAMAF"
    threads : 1
    resources:
        mem_mb = 10000
    shell:
        'oncotator --input_format=VCF --output_format=TCGAMAF --tx-mode EFFECT --db-dir={input.ONCOTATOR_DB} {input.interval_vcf} {output.MAF} hg19 2> {log}'
        #'set +u; source {params.activate_oncotator}; set -u; oncotator --input_format=VCF --output_format=TCGAMAF --tx-mode EFFECT --db-dir={input.ONCOTATOR_DB} {input.interval_vcf} {output.MAF} hg19; deactivate'

# concatenate oncotator T_only
rule concatenate_oncotator_tumor_only:
    input:
        maf = expand("oncotator_T_tmp/{{tsample}}_tumor_only_ON_{mutect_interval}_annotated_T.TCGAMAF", mutect_interval=mutect_intervals)
    output:
        concatened_oncotator = "oncotator_T/{tsample}_tumor_only_annotated_T.TCGAMAF",
        tmp_list = temp("oncotator_T_tmp/{tsample}_T_oncotator_tmp.list")
    params:
        queue = "shortq",
        merge_oncotator = config["APP_MERGE_ONCOTATOR"]
    threads : 1
    resources:
        mem_mb = 10000
    log:
        "logs/merge_oncotator/{tsample}_tumor_only_annotated_T.log"
    shell :
        "ls -1a oncotator_T_tmp/{wildcards.tsample}_tumor_only_ON_*_annotated_T.TCGAMAF > oncotator_T_tmp/{wildcards.tsample}_T_oncotator_tmp.list && "
        "python2.7  {params.merge_oncotator} {output.tmp_list} {output.concatened_oncotator} 2> {log}" 

## A rule to simplify oncotator output on tumor only samples
rule oncotator_reformat_tumor_only:
    input:
        maf="oncotator_T/{tsample}_tumor_only_annotated_T.TCGAMAF"
    output:
        maf ="oncotator_T_maf/{tsample}_tumor_only_T_selection.TCGAMAF",
        tsv ="oncotator_T_tsv/{tsample}_tumor_only_T.tsv"
    log:
        "logs/oncotator/{tsample}_tumor_only_T_selection.txt"
    params:
        queue = "shortq",
        oncotator_extract_Tonly = config["APP_ONCOTATOR_EXTRACT_TUMOR_ONLY"]
    threads : 1
    resources:
        mem_mb = 10000
    shell:
        'python2.7 {params.oncotator_extract_Tonly} {input.maf} {output.maf} {output.tsv} 2> {log}'

## A rule to simplify oncotator output on tumor only samples
rule oncotator_with_pileup_tumor_only:
    input:
        tsv = "oncotator_T_tsv/{tsample}_tumor_only_T.tsv",
        pileup = "pileup_T/{tsample}_tumor_only_T.pileup.gz"
    output:
        tsv = "oncotator_T_tsv_pileup/{tsample}_tumor_only_T_with_pileup.tsv"
    log:
        "logs/oncotator/{tsample}_tumor_only_T_selection_with_pileup.txt"
    params:
        queue = "shortq",
        oncotator_cross_pileup = config["APP_ONCOTATOR_X_PILEUP"]
    threads : 1
    resources:
        mem_mb = 500
    shell:
        'python {params.oncotator_cross_pileup} {input.pileup} {input.tsv} {output.tsv}'

## A rule to simplify oncotator output on tumor only samples
rule oncotator_with_COSMIC_tumor_only:
    input:
        tsv = "oncotator_T_tsv_pileup/{tsample}_tumor_only_T_with_pileup.tsv"
    output:
        tsv = "oncotator_T_tsv_COSMIC/{tsample}_tumor_only_T_with_COSMIC.tsv"
    log:
        "logs/oncotator/{tsample}_tumor_only_T_selection_with_COSMIC.txt"
    params:
        queue = "shortq",
        oncotator_cross_cosmic = config["APP_ONCOTATOR_X_COSMIC_T_ONLY"],
        cosmic_mutation = config["COSMIC_MUTATION"],
        cancer_census_oncogene = config["CANCER_CENSUS_ONCOGENE"],
        cancer_census_tumorsupressor = config["CANCER_CENSUS_TUMORSUPRESSOR"]
    threads : 1
    resources:
        mem_mb = 10000
    shell:
        'python2.7 {params.oncotator_cross_cosmic}  {input.tsv} {output.tsv} {params.cosmic_mutation} {params.cancer_census_oncogene} {params.cancer_census_tumorsupressor} 2> {log}'

