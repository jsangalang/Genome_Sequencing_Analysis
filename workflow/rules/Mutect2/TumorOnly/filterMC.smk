## A rule to filter variant call, from mutect tumor only
rule filter_mutect_calls_tumor_only:
    input :
        Mutect2_vcf = "Mutect2_T/{tsample}_tumor_only_T.vcf.gz",
        Mutect2_stats = "Mutect2_T/{tsample}_tumor_only_T.vcf.gz.stats",
        contamination_table = "cross_sample_contamination/{tsample}_calculatecontamination.table",
        MUTECT_FILTER_REF = config["MUTECT_FILTER_REF"],
        index = config["GENOME_FASTA"]
    output:
        VCF = "Mutect2_T/{tsample}_tumor_only_filtered_T.vcf.gz",
        INDEX = "Mutect2_T/{tsample}_tumor_only_filtered_T.vcf.gz.tbi"
    params:
        queue = "mediumq",
        gatk = config["APP_GATK"]
    log:
        "logs/filter_Mutect2_T/{tsample}_tumor_only_filtered_T.vcf.gz.log"
    threads : 1
    resources:
        mem_mb = 100000
    shell:
        "{params.gatk} --java-options \"-Xmx40g  -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp \" FilterMutectCalls"
        " -V {input.Mutect2_vcf}"
        " -R {input.index}"
        " --contamination-table {input.contamination_table}"
        " -O {output.VCF} 2> {log}"
        
