# A rule to generate a bed from mutect2 vcf, on tumor only with panel of normal
rule get_variant_bed_tumor_only_pon:
    input:
        Mutect2_vcf = "Mutect2_Tp/{tsample}_PON_{panel_of_normal}_twicefiltered_Tp.vcf.gz",
        Mutect2_vcf_index = "Mutect2_Tp/{tsample}_PON_{panel_of_normal}_twicefiltered_Tp.vcf.gz.tbi",
    output:
        BED = "variant_bed_Tp/{tsample}_PON_{panel_of_normal}_Tp.bed"
    log:
        "logs/variant_bed_Tp/{tsample}_PON_{panel_of_normal}_Tp.bed.txt"
    params:
        queue = "mediumq",
        vcf2bed = config["APP_VCF2BED"] 
    threads : 1
    resources:
        mem_mb = 5000
    shell:
        'zcat {input.Mutect2_vcf} | python2 {params.vcf2bed} - > {output.BED} 2> {log}'

# Run samtools mpileup, on tumor only with panel of normal
rule samtools_mpileup_tumor_only_pon:
    input:
        BED = "variant_bed_Tp/{tsample}_PON_{panel_of_normal}_Tp.bed",
        GENOME_REF_FASTA = config["GENOME_FASTA"],
        BAM = "bam/{tsample}.nodup.recal.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{tsample}.recal.bam",
        BAI = "bam/{tsample}.nodup.recal.bam.bai" if config["REMOVE_DUPLICATES"]=="True" else "bam/{tsample}.recal.bam.bai"
    output:
        PILEUP = "pileup_Tp/{tsample}_PON_{panel_of_normal}_Tp.pileup.gz"
    log:
        "logs/pileup_Tp/{tsample}_PON_{panel_of_normal}_Tp.pileup.txt"
    params:
        queue = "mediumq",
        samtools = config["APP_SAMTOOLS"]
    threads : 1
    resources:
        mem_mb = 5000
    shell:
        '{params.samtools} mpileup -a -B -l {input.BED} -f {input.GENOME_REF_FASTA} {input.BAM} | gzip - > {output.PILEUP} 2> {log}'
        
## A rule to split mutect2 results in pieces 
rule split_Mutect2_tumor_only_pon:
    input:
        Mutect2_vcf = "Mutect2_Tp/{tsample}_PON_{panel_of_normal}_twicefiltered_Tp.vcf.gz",
        vcf_index = "Mutect2_Tp/{tsample}_PON_{panel_of_normal}_twicefiltered_Tp.vcf.gz.tbi",
        interval = config["MUTECT_INTERVAL_DIR"] + "/{interval}.bed"
    output:
        interval_vcf_bcftools = temp("Mutect2_Tp_oncotator_tmp/{tsample}_PON_{panel_of_normal}_Tp_ON_{interval}_bcftools.vcf.gz"),
        interval_vcf = temp("Mutect2_Tp_oncotator_tmp/{tsample}_PON_{panel_of_normal}_Tp_ON_{interval}.vcf.gz")
    params:
        queue = "shortq",
        bcftools = config["APP_BCFTOOLS"],
        reformat_script = config["APP_REFORMAT_MUTECT2"]
    log:
        "logs/Mutect2_Tp_oncotator_tmp/{tsample}_PON_{panel_of_normal}_Tp_ON_{interval}.vcf.log"
    threads : 1
    resources:
        mem_mb = 2000
    shell:
        '{params.bcftools} view -l 9 -R {input.interval} -o {output.interval_vcf_bcftools} {input.Mutect2_vcf} 2> {log}  &&'
        ' python {params.reformat_script} {output.interval_vcf_bcftools} {output.interval_vcf} 2>> {log}'        
        
# A rule to annotate mutect2 tumor only with panel of normal results with oncotator 
rule oncotator_tumor_only_pon:
    input:
        ONCOTATOR_DB = config["ONCOTATOR_DB"],
        interval_vcf = "Mutect2_Tp_oncotator_tmp/{tsample}_PON_{panel_of_normal}_Tp_ON_{interval}.vcf.gz"
    output:
        MAF = temp("oncotator_Tp_tmp/{tsample}_PON_{panel_of_normal}_ON_{interval}_annotated_Tp.TCGAMAF")
    params:
        queue = "mediumq",
        #activate_oncotator = config["ONCOTATOR_ENV"],
        oncotator = config["APP_ONCOTATOR"]
    log:
        "logs/oncotator_Tp_tmp/{tsample}_PON_{panel_of_normal}_ON_{interval}_annotated_Tp.TCGAMAF"
    threads : 1
    resources:
        mem_mb = 10000
    shell:
        'oncotator --input_format=VCF --output_format=TCGAMAF --tx-mode EFFECT --db-dir={input.ONCOTATOR_DB} {input.interval_vcf} {output.MAF} hg19 2> {log}'
        #'set +u; source {params.activate_oncotator}; set -u; oncotator --input_format=VCF --output_format=TCGAMAF --tx-mode EFFECT --db-dir={input.ONCOTATOR_DB} {input.interval_vcf} {output.MAF} hg19; deactivate'

# concatenate oncotator T_only_pon
rule concatenate_oncotator_tumor_only_pon:
    input:
        maf = expand("oncotator_Tp_tmp/{{tsample}}_PON_{{panel_of_normal}}_ON_{mutect_interval}_annotated_Tp.TCGAMAF", mutect_interval=mutect_intervals)
    output:
        concatened_oncotator = "oncotator_Tp/{tsample}_PON_{panel_of_normal}_annotated_Tp.TCGAMAF",
        tmp_list = temp("oncotator_Tp_tmp/{tsample}_PON_{panel_of_normal}_Tp_oncotator_tmp.list")
    params:
        queue = "shortq",
        merge_oncotator = config["APP_MERGE_ONCOTATOR"]
    threads : 1
    resources:
        mem_mb = 10000
    log:
        "logs/merge_oncotator/{tsample}_PON_{panel_of_normal}_annotated_Tp.log"
    shell :
        "ls -1a oncotator_Tp_tmp/{wildcards.tsample}_PON_{wildcards.panel_of_normal}_ON_*_annotated_Tp.TCGAMAF > oncotator_Tp_tmp/{wildcards.tsample}_PON_{wildcards.panel_of_normal}_Tp_oncotator_tmp.list && "
        "python2.7  {params.merge_oncotator} {output.tmp_list} {output.concatened_oncotator} 2> {log}" 

## A rule to simplify oncotator output on tumor only samples with pon
rule oncotator_reformat_tumor_only_pon:
    input:
        maf="oncotator_Tp/{tsample}_PON_{panel_of_normal}_annotated_Tp.TCGAMAF"
    output:
        maf ="oncotator_Tp_maf/{tsample}_PON_{panel_of_normal}_Tp_selection.TCGAMAF",
        tsv ="oncotator_Tp_tsv/{tsample}_PON_{panel_of_normal}_Tp.tsv"
    log:
        "logs/oncotator/{tsample}_PON_{panel_of_normal}_Tp_selection.txt"
    params:
        queue = "shortq",
        oncotator_extract_Tonly = config["APP_ONCOTATOR_EXTRACT_TUMOR_ONLY"]
    threads : 1
    resources:
        mem_mb = 10000
    shell:
        'python2.7 {params.oncotator_extract_Tonly} {input.maf} {output.maf} {output.tsv} 2> {log}'

## A rule to simplify oncotator output on tumor only samples with pon
rule oncotator_with_pileup_tumor_only_pon:
    input:
        tsv = "oncotator_Tp_tsv/{tsample}_PON_{panel_of_normal}_Tp.tsv",
        pileup = "pileup_Tp/{tsample}_PON_{panel_of_normal}_Tp.pileup.gz"
    output:
        tsv = "oncotator_Tp_tsv_pileup/{tsample}_PON_{panel_of_normal}_Tp_with_pileup.tsv"
    log:
        "logs/oncotator/{tsample}_PON_{panel_of_normal}_Tp_with_pileup.txt"
    params:
        queue = "shortq",
        oncotator_cross_pileup = config["APP_ONCOTATOR_X_PILEUP"]
    threads : 1
    resources:
        mem_mb = 500
    shell:
        'python {params.oncotator_cross_pileup} {input.pileup} {input.tsv} {output.tsv}'
        
## A rule to simplify oncotator output on tumor only samples with pon
rule oncotator_with_COSMIC_tumor_only_pon:
    input:
        tsv = "oncotator_Tp_tsv_pileup/{tsample}_PON_{panel_of_normal}_Tp_with_pileup.tsv"
    output:
        tsv = "oncotator_Tp_tsv_COSMIC/{tsample}_PON_{panel_of_normal}_Tp_with_COSMIC.tsv"
    log:
        "logs/oncotator/{tsample}_PON_{panel_of_normal}_Tp_with_COSMIC.txt"
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
