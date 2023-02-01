# A rule to generate a bed from mutect2 vcf, on tumor only 
rule extract_exom_mutect2_tumor_only:
    input:
        Mutect2_vcf = "Mutect2_T/{tsample}_tumor_only_twicefiltered_T.vcf.gz",
        Mutect2_vcf_index = "Mutect2_T/{tsample}_tumor_only_twicefiltered_T.vcf.gz.tbi",
        exom_bed = config["EXOM_BED"]
    output:
        exom_Mutect2 = temp("Mutect2_T_exom/{tsample}_tumor_only_twicefiltered_T_exom_unsorted.vcf.gz")
    log:
        "logs/Mutect2_T_exom/{tsample}_tumor_only_T.vcf.txt"
    params:
        queue = "mediumq",
        bcftools = config["APP_BCFTOOLS"] 
    threads : 1
    resources:
        mem_mb = 5000
    shell:
        '{params.bcftools} view -l 9 -R {input.exom_bed} -o {output.exom_Mutect2} {input.Mutect2_vcf} 2> {log}'
 
# A rule to sort exom vcf
rule sort_exom_mutect2_tumor_only:
    input:
        Mutect2_vcf = "Mutect2_T_exom/{tsample}_tumor_only_twicefiltered_T_exom_unsorted.vcf.gz"
    output:
        exom_Mutect2 = "Mutect2_T_exom/{tsample}_tumor_only_twicefiltered_T_exom.vcf.gz"
    log:
        "logs/Mutect2_T_exom/{tsample}_tumor_only_T_sort.txt"
    params:
        queue = "shortq",
        vcfSort = config["APP_VCFSORT"]
    threads : 1
    resources:
        mem_mb = 5000
    shell:
        'bgzip -d {input.Mutect2_vcf} 2>> log && '
        '{params.vcfSort} Mutect2_T_exom/{wildcards.tsample}_tumor_only_twicefiltered_T_exom_unsorted.vcf > Mutect2_T_exom/{wildcards.tsample}_tumor_only_twicefiltered_T_exom.vcf  2>> log && '
        'bgzip Mutect2_T_exom/{wildcards.tsample}_tumor_only_twicefiltered_T_exom.vcf 2>> log'
        
# A rule to generate a bed from mutect2 vcf, on tumor only 
rule index_exom_mutect2_tumor_only:
    input:
        exom_Mutect2 = "Mutect2_T_exom/{tsample}_tumor_only_twicefiltered_T_exom.vcf.gz"
    output:
        exom_Mutect2 = "Mutect2_T_exom/{tsample}_tumor_only_twicefiltered_T_exom.vcf.gz.tbi"
    log:
        "logs/Mutect2_T_exom/{tsample}_tumor_only_T_index.txt"
    params:
        queue = "shortq",
        gatk = config["APP_GATK"]
    threads : 1
    conda: "pipeline_GATK_2.1.4_V2"
    resources:
        mem_mb = 1000
    shell:
        'gatk IndexFeatureFile -F {input.exom_Mutect2} 2> {log}'
        
# A rule to generate a bed from mutect2 vcf, on tumor only 
rule get_variant_bed_tumor_only_exom:
    input:
        Mutect2_vcf = "Mutect2_T_exom/{tsample}_tumor_only_twicefiltered_T_exom.vcf.gz",
        Mutect2_vcf_index = "Mutect2_T_exom/{tsample}_tumor_only_twicefiltered_T_exom.vcf.gz.tbi"
    output:
        BED = "variant_bed_T_exom/{tsample}_tumor_only_T_exom.bed"
    log:
        "logs/variant_bed_T_exom/{tsample}_tumor_only_T_exom.bed.txt"
    params:
        queue = "mediumq",
        vcf2bed = config["APP_VCF2BED"] 
    threads : 1
    resources:
        mem_mb = 5000
    shell:
        'zcat {input.Mutect2_vcf} | python2 {params.vcf2bed} - > {output.BED} 2> {log}'

# Run samtools mpileup, on tumor only 
rule samtools_mpileup_tumor_only_exom:
    input:
        BED = "variant_bed_T_exom/{tsample}_tumor_only_T_exom.bed",
        GENOME_REF_FASTA = config["GENOME_FASTA"],
        BAM = "bam/{tsample}.nodup.recal.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{tsample}.recal.bam",
        BAI = "bam/{tsample}.nodup.recal.bam.bai" if config["REMOVE_DUPLICATES"]=="True" else "bam/{tsample}.recal.bam.bai"
    output:
        PILEUP = "pileup_T_exom/{tsample}_tumor_only_T_exom.pileup.gz"
    log:
        "logs/pileup_T_exom/{tsample}_tumor_only_T_exom.pileup.txt"
    params:
        queue = "mediumq",
        samtools = config["APP_SAMTOOLS"]
    threads : 1
    resources:
        mem_mb = 5000
    shell:
        '{params.samtools} mpileup -a -B -l {input.BED} -f {input.GENOME_REF_FASTA} {input.BAM} | gzip - > {output.PILEUP} 2> {log}'

# A rule to annotate mutect2 tumor only results with oncotator 
rule oncotator_tumor_only_exom:
    input:
        ONCOTATOR_DB = config["ONCOTATOR_DB"],
        Mutect2_vcf = "Mutect2_T_exom/{tsample}_tumor_only_twicefiltered_T_exom.vcf.gz",
        Mutect2_vcf_index = "Mutect2_T_exom/{tsample}_tumor_only_twicefiltered_T_exom.vcf.gz.tbi"
    output:
        MAF="oncotator_T_exom/{tsample}_tumor_only_annotated_T_exom.TCGAMAF"
    params:
        queue = "mediumq",
        #activate_oncotator = config["ONCOTATOR_ENV"],
        oncotator = config["APP_ONCOTATOR"]
    log:
        "logs/oncotator_T_exom/{tsample}_tumor_only_annotated_T_exom.TCGAMAF"
    threads : 1
    resources:
        mem_mb = 100000
    shell:
        '{params.oncotator} --input_format=VCF --output_format=TCGAMAF --tx-mode EFFECT --db-dir={input.ONCOTATOR_DB} {input.Mutect2_vcf} {output.MAF} hg19 2> {log}'
        #'set +u; source {params.activate_oncotator}; set -u; oncotator --input_format=VCF --output_format=TCGAMAF --tx-mode EFFECT --db-dir={input.ONCOTATOR_DB} {input.Mutect2_vcf} {output.MAF} hg19; deactivate'

## A rule to simplify oncotator output on tumor only samples
rule oncotator_reformat_tumor_only_exom:
    input:
        maf="oncotator_T_exom/{tsample}_tumor_only_annotated_T_exom.TCGAMAF"
    output:
        maf ="oncotator_T_maf_exom/{tsample}_tumor_only_T_selection_exom.TCGAMAF",
        tsv ="oncotator_T_tsv_exom/{tsample}_tumor_only_T_exom.tsv"
    log:
        "logs/oncotator_exom/{tsample}_tumor_only_T_selection_exom.txt"
    params:
        queue = "shortq",
        oncotator_extract_Tonly = config["APP_ONCOTATOR_EXTRACT_TUMOR_ONLY"]
    threads : 1
    resources:
        mem_mb = 10000
    shell:
        'python2.7 {params.oncotator_extract_Tonly} {input.maf} {output.maf} {output.tsv} 2> {log}'

## A rule to simplify oncotator output on tumor only samples
rule oncotator_with_pileup_tumor_only_exom:
    input:
        tsv = "oncotator_T_tsv_exom/{tsample}_tumor_only_T_exom.tsv",
        pileup = "pileup_T_exom/{tsample}_tumor_only_T_exom.pileup.gz"
    output:
        tsv = "oncotator_T_tsv_pileup_exom/{tsample}_tumor_only_T_with_pileup_exom.tsv"
    log:
        "logs/oncotator_exom/{tsample}_tumor_only_T_selection_with_pileup_exom.txt"
    params:
        queue = "shortq",
        oncotator_cross_pileup = config["APP_ONCOTATOR_X_PILEUP"]
    threads : 1
    resources:
        mem_mb = 500
    shell:
        'python {params.oncotator_cross_pileup} {input.pileup} {input.tsv} {output.tsv}'

## A rule to simplify oncotator output on tumor only samples
rule oncotator_with_COSMIC_tumor_only_exom:
    input:
        tsv = "oncotator_T_tsv_pileup_exom/{tsample}_tumor_only_T_with_pileup_exom.tsv"
    output:
        tsv = "oncotator_T_tsv_COSMIC_exom/{tsample}_tumor_only_T_with_COSMIC_exom.tsv"
    log:
        "logs/oncotator_exom/{tsample}_tumor_only_T_selection_with_COSMIC_exom.txt"
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

