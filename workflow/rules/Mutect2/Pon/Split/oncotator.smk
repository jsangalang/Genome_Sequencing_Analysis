# A rule to generate a bed from mutect2 vcf, on tumor versus normal with panel of normals
rule get_variant_bed_pon:
    input:
        Mutect2_vcf = "Mutect2_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_twicefiltered_TvNp.vcf.gz"
    output:
        BED = "variant_bed_TvN/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvN.bed"
    log:
        "logs/variant_bed_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvN.bed.txt"
    params:
        queue = "mediumq",
        vcf2bed = config["vcf2bed"]["app"]
    threads : 1
    resources:
        mem_mb = 10240
    shell:
        'zcat {input.Mutect2_vcf} | python2 {params.vcf2bed} - > {output.BED} 2> {log}'

# Run samtools mpileup, on tumor versus normal with panel of normals
rule samtools_mpileup_pon:
    input:
        BED = "variant_bed_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp.bed",
        BAM = "bam/{tsample}.nodup.recal.bam" if config["REMOVE_DUPLICATES"] == True else "bam/{tsample}.recal.bam",
        BAI = "bam/{tsample}.nodup.recal.bam.bai" if config["REMOVE_DUPLICATES"] == True else "bam/{tsample}.recal.bam.bai"
    output:
        PILEUP = "pileup_TvN/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvN.pileup.gz"
    log:
        "logs/pileup_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvN.pileup.txt"
    params:
        queue = "mediumq",
        samtools = config["samtools"]["app"],
        GENOME_REF_FASTA = config["gatk"][config["samples"]]["genome_fasta"],
    threads : 1
    resources:
        mem_mb = 10240
    shell:
        '{params.samtools} mpileup -a -B -l {input.BED} -f {params.GENOME_REF_FASTA} {input.BAM} | gzip - > {output.PILEUP} 2> {log}'

## A rule to split mutect2 results in pieces 
rule split_Mutect2_pon:
    input:
        Mutect2_vcf = "Mutect2_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_twicefiltered_TvNp.vcf.gz",
        vcf_index = "Mutect2_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_twicefiltered_TvNp.vcf.gz.tbi",
    output:
        interval_vcf_bcftools = temp("Mutect2_TvNp_oncotator_tmp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp_ON_{interval}_bcftools.vcf.gz"),
        interval_vcf = temp("Mutect2_TvNp_oncotator_tmp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp_ON_{interval}.vcf.gz")
    params:
        queue = "shortq",
        bcftools = config["bcftools"]["app"],
        reformat = config["gatk"]["scripts"]["reformat_mutect2"],
        interval = config["gatk"][config["samples"]][config["seq_type"]]["mutect_interval_dir"] + "/{interval}.bed"
    log:
        "logs/Mutect2_TvNp_oncotator_tmp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp_ON_{interval}.vcf.log"
    threads : 1
    resources:
        mem_mb = 20480
    shell:
        '{params.bcftools} view -l 9 -R {input.interval} -o {output.interval_vcf_bcftools} {input.Mutect2_vcf} 2> {log} &&'
        ' python {params.reformat_script} {output.interval_vcf_bcftools} {output.interval_vcf} 2>> {log}'

# A rule to annotate mutect2 tumor versus normal and panel of normal results with oncotator  
rule oncotator_pon:
    input:
        interval_vcf = "Mutect2_TvNp_oncotator_tmp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp_ON_{interval}.vcf.gz"
    output:
        MAF = temp("oncotator_TvNp_tmp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_ON_{interval}_annotated_TvNp.TCGAMAF")
    params:
        queue = "mediumq",
        oncotator = config["oncotator"]["app"],
        oncotator_db = config["oncotator"][config["samples"]]["DB"],
    log:
        "logs/oncotator_TvNp_tmp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_ON_{interval}_annotated_TvNp.TCGAMAF.log"
    threads : 1
    resources:
        mem_mb = 10240
    shell:
        '{params.oncotator} --input_format=VCF --output_format=TCGAMAF --tx-mode EFFECT --db-dir={params.oncotator_db} {input.interval_vcf} {output.MAF} hg19 2> {log}'

# concatenate oncotator TvN_pon
rule concatenate_oncotator_pon:
    input:
        maf = expand("oncotator_TvNp_tmp/{{tsample}}_Vs_{{nsample}}_PON_{{panel_of_normal}}_ON_{mutect_interval}_annotated_TvNp.TCGAMAF", mutect_interval=mutect_intervals)
    output:
        concatened_oncotator = "oncotator_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_annotated_TvNp.TCGAMAF",
        tmp_list = temp("oncotator_TvNp_tmp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp_oncotator_tmp.list")
    params:
        queue = "shortq",
        merge_oncotator = config["APP_MERGE_ONCOTATOR"]
    threads : 1
    resources:
        mem_mb = 10000
    log:
        "logs/merge_oncotator/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp.vcf.log"
    shell :
        "ls -1a oncotator_TvNp_tmp/{wildcards.tsample}_Vs_{wildcards.nsample}_PON_{wildcards.panel_of_normal}_ON_*_annotated_TvNp.TCGAMAF > oncotator_TvNp_tmp/{wildcards.tsample}_Vs_{wildcards.nsample}_PON_{wildcards.panel_of_normal}_TvNp_oncotator_tmp.list && "
        "python2.7  {params.merge_oncotator} {output.tmp_list} {output.concatened_oncotator} 2> {log}" 

## A rule to simplify oncotator output on tumor vs normal samples with panel of normal
rule oncotator_reformat_TvN_pon:
    input:
        maf="oncotator_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_annotated_TvNp.TCGAMAF"
    output:
        maf ="oncotator_TvNp_maf/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp_selection.TCGAMAF",
        tsv ="oncotator_TvNp_tsv/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp.tsv"
    log:
        "logs/oncotator_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_annotated_TvNp_selection.txt"
    params:
        queue = "shortq",
        oncotator_extract_TvN = config["APP_ONCOTATOR_EXTRACT_TUMOR_VS_NORMAL"]
    threads : 1
    resources:
        mem_mb = 10000
    shell:
        'python2.7 {params.oncotator_extract_TvN} {input.maf} {output.maf} {output.tsv} 2> {log}'

## A rule to simplify oncotator output on tumor vs normal samples with panel of normal
rule oncotator_with_pileup_TvN_pon:
    input:
        tsv = "oncotator_TvNp_tsv/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp.tsv",
        pileup = "pileup_TvN/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvN.pileup.gz"
    output:
        tsv = "oncotator_TvNp_tsv_pileup/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp_with_pileup.tsv"
    log:
        "logs/oncotator/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_annotated_TvNp_with_pileup.txt"
    params:
        queue = "shortq",
        oncotator_cross_pileup = config["APP_ONCOTATOR_X_PILEUP"]
    threads : 1
    resources:
        mem_mb = 500
    shell:
        'python {params.oncotator_cross_pileup} {input.pileup} {input.tsv} {output.tsv}'

## A rule to simplify oncotator output on tumor vs normal samples with panel of normal
rule oncotator_with_COSMIC_TvN_pon:
    input:
        tsv = "oncotator_TvNp_tsv_pileup/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp_with_pileup.tsv"
    output:
        tsv = "oncotator_TvNp_tsv_COSMIC/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp_with_COSMIC.tsv"
    log:
        "logs/oncotator/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_annotated_TvNp_with_COSMIC.txt"
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

