# A rule to generate a bed from mutect2 vcf  
rule get_variant_bed:
    input:
        Mutect2_vcf = "Mutect2_TvN/{tsample}_Vs_{nsample}_twicefiltered_TvN.vcf.gz"
    output:
        BED = "variant_bed_TvN/{tsample}_Vs_{nsample}_TvN.bed"
    log:
        "logs/variant_bed_TvN/{tsample}_Vs_{nsample}_TvN.bed.txt"
    params:
        queue = "mediumq",
        vcf2bed = config["APP_VCF2BED"]        
    threads : 1
    resources:
        mem_mb = 5000
    shell:
        'zcat {input.Mutect2_vcf} | python2 {params.vcf2bed} - > {output.BED} 2> {log}'

## Run samtools mpileup 
rule samtools_mpileup:
    input:
        BED = "variant_bed_TvN/{tsample}_Vs_{nsample}_TvN.bed",
        GENOME_REF_FASTA = config["GENOME_FASTA"],
        BAM = "bam/{tsample}.nodup.recal.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{tsample}.recal.bam",
        BAI = "bam/{tsample}.nodup.recal.bam.bai" if config["REMOVE_DUPLICATES"]=="True" else "bam/{tsample}.recal.bam.bai"
    output:
        PILEUP = "pileup_TvN/{tsample}_Vs_{nsample}_TvN.pileup.gz"
    log:
        "logs/pileup_TvN/{tsample}_Vs_{nsample}_TvN.pileup.txt"
    params:
        queue = "mediumq",
        samtools = config["APP_SAMTOOLS"]
    threads : 1
    resources:
        mem_mb = 5000
    shell:
        '{params.samtools} mpileup -a -B -l {input.BED} -f {input.GENOME_REF_FASTA} {input.BAM} | gzip - > {output.PILEUP} 2> {log}'

## A rule to split mutect2 results in pieces 
rule split_Mutect2:
    input:
        Mutect2_vcf = "Mutect2_TvN/{tsample}_Vs_{nsample}_twicefiltered_TvN.vcf.gz",
        vcf_index = "Mutect2_TvN/{tsample}_Vs_{nsample}_twicefiltered_TvN.vcf.gz.tbi",
        interval = config["MUTECT_INTERVAL_DIR"] + "/{interval}.bed"
    output:
        interval_vcf_bcftools = temp("Mutect2_TvN_oncotator_tmp/{tsample}_Vs_{nsample}_TvN_ON_{interval}_bcftools.vcf.gz"),
        interval_vcf = temp("Mutect2_TvN_oncotator_tmp/{tsample}_Vs_{nsample}_TvN_ON_{interval}.vcf.gz")
    params:
        queue = "shortq",
        bcftools = config["APP_BCFTOOLS"],
        reformat_script = config["APP_REFORMAT_MUTECT2"]
    log:
        "logs/Mutect2_TvN_oncotator_tmp/{tsample}_Vs_{nsample}_TvN_ON_{interval}.vcf.log"
    threads : 1
    resources:
        mem_mb = 2000
    shell: 
        '{params.bcftools} view -l 9 -R {input.interval} -o {output.interval_vcf_bcftools} {input.Mutect2_vcf} 2> {log} &&'
        ' python {params.reformat_script} {output.interval_vcf_bcftools} {output.interval_vcf} 2>> {log}'

# A rule to annotate mutect2 tumor versus normal results with oncotator  
rule oncotator:
    input:
        ONCOTATOR_DB = config["ONCOTATOR_DB"],
        interval_vcf = "Mutect2_TvN_oncotator_tmp/{tsample}_Vs_{nsample}_TvN_ON_{interval}.vcf.gz"
    output:
        MAF = temp("oncotator_TvN_tmp/{tsample}_Vs_{nsample}_ON_{interval}_annotated_TvN.TCGAMAF")
    params:
        queue = "mediumq",
        #activate_oncotator = config["ONCOTATOR_ENV"],
        oncotator = config["APP_ONCOTATOR"]
    log:
        "logs/oncotator_TvN_tmp/{tsample}_Vs_{nsample}_ON_{interval}_annotated_TvN.TCGAMAF"
    threads : 1
    resources:
        mem_mb = 10000
    shell:
        'oncotator --input_format=VCF --output_format=TCGAMAF --tx-mode EFFECT --db-dir={input.ONCOTATOR_DB} {input.interval_vcf} {output.MAF} hg19 2> {log}'
        #'set +u; source {params.activate_oncotator}; set -u; oncotator --input_format=VCF --output_format=TCGAMAF --tx-mode EFFECT --db-dir={input.ONCOTATOR_DB} {input.interval_vcf} {output.MAF} hg19; deactivate;'

# concatenate oncotator TvN
rule concatenate_oncotator:
    input:
        maf = expand("oncotator_TvN_tmp/{{tsample}}_Vs_{{nsample}}_ON_{mutect_interval}_annotated_TvN.TCGAMAF", mutect_interval=mutect_intervals)
    output:
        concatened_oncotator = "oncotator_TvN/{tsample}_Vs_{nsample}_annotated_TvN.TCGAMAF",
        tmp_list = temp( "oncotator_TvN_tmp/{tsample}_Vs_{nsample}_TvN_oncotator_tmp.list")
    params:
        queue = "shortq",
        merge_oncotator = config["APP_MERGE_ONCOTATOR"]
    threads : 1
    resources:
        mem_mb = 10000
    log:
        "logs/merge_oncotator/{tsample}_Vs_{nsample}_TvN.vcf.log"
    shell :
        "ls -1a oncotator_TvN_tmp/{wildcards.tsample}_Vs_{wildcards.nsample}_ON_*_annotated_TvN.TCGAMAF > oncotator_TvN_tmp/{wildcards.tsample}_Vs_{wildcards.nsample}_TvN_oncotator_tmp.list && "
        "python2.7  {params.merge_oncotator} {output.tmp_list} {output.concatened_oncotator} 2> {log}"
        
## A rule to simplify oncotator output on tumor vs normal samples
rule oncotator_reformat_TvN:
    input:
        maf="oncotator_TvN/{tsample}_Vs_{nsample}_annotated_TvN.TCGAMAF"
    output:
        maf ="oncotator_TvN_maf/{tsample}_Vs_{nsample}_TvN_selection.TCGAMAF",
        tsv ="oncotator_TvN_tsv/{tsample}_Vs_{nsample}_TvN.tsv"
    log:
        "logs/oncotator/{tsample}_Vs_{nsample}_TvN_selection.txt"
    params:
        queue = "shortq",
        oncotator_extract_TvN = config["APP_ONCOTATOR_EXTRACT_TUMOR_VS_NORMAL"]
    threads : 1
    resources:
        mem_mb = 10000
    shell:
        'python2.7 {params.oncotator_extract_TvN} {input.maf} {output.maf} {output.tsv} 2> {log}'

## A rule to cross oncotator output on tumor vs normal samples with pileup information
rule oncotator_with_pileup_TvN:
    input:
        tsv = "oncotator_TvN_tsv/{tsample}_Vs_{nsample}_TvN.tsv",
        pileup = "pileup_TvN/{tsample}_Vs_{nsample}_TvN.pileup.gz"
    output:
        tsv = "oncotator_TvN_tsv_pileup/{tsample}_Vs_{nsample}_TvN_with_pileup.tsv"
    log:
        "logs/oncotator/{tsample}_Vs_{nsample}_TvN_with_pileup.txt"
    params:
        queue = "shortq",
        oncotator_cross_pileup = config["APP_ONCOTATOR_X_PILEUP"]
    threads : 1
    resources:
        mem_mb = 500
    shell:
        'python {params.oncotator_cross_pileup} {input.pileup} {input.tsv} {output.tsv}'

## A rule to cross oncotator output on tumor vs normal samples with COSMIC information
rule oncotator_with_COSMIC_TvN:
    input:
        tsv = "oncotator_TvN_tsv_pileup/{tsample}_Vs_{nsample}_TvN_with_pileup.tsv"
    output:
        tsv = "oncotator_TvN_tsv_COSMIC/{tsample}_Vs_{nsample}_TvN_with_COSMIC.tsv"
    log:
        "logs/oncotator/{tsample}_Vs_{nsample}_TvN_with_COSMIC.txt"
    params:
        queue = "shortq",
        oncotator_cross_cosmic = config["APP_ONCOTATOR_X_COSMIC_T_N"],
        cosmic_mutation = config["COSMIC_MUTATION"],
        cancer_census_oncogene = config["CANCER_CENSUS_ONCOGENE"],
        cancer_census_tumorsupressor = config["CANCER_CENSUS_TUMORSUPRESSOR"]
    threads : 1
    resources:
        mem_mb = 10000
    shell:
        'python2.7 {params.oncotator_cross_cosmic} {input.tsv} {output.tsv} {params.cosmic_mutation} {params.cancer_census_oncogene} {params.cancer_census_tumorsupressor} 2> {log}'

