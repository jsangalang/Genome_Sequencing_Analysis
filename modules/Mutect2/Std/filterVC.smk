## A rule to filter variant call, from mutect tumor Vs normal
rule filter_mutect_calls:
    input :
        Mutect2_vcf = "Mutect2_TvN/{tsample}_Vs_{nsample}_TvN.vcf.gz",
        Mutect2_stats = "Mutect2_TvN/{tsample}_Vs_{nsample}_TvN.vcf.gz.stats",
        contamination_table = "cross_sample_contamination/{tsample}_calculatecontamination.table",
        index = config["GENOME_FASTA"]
    output:
        VCF = "Mutect2_TvN/{tsample}_Vs_{nsample}_filtered_TvN.vcf.gz",
        INDEX = "Mutect2_TvN/{tsample}_Vs_{nsample}_filtered_TvN.vcf.gz.tbi"
    params:
        queue = "mediumq",
        gatk = config["APP_GATK"]
    log:
        "logs/filter_Mutect2_TvN/{tsample}_Vs_{nsample}_filtered_TvN.vcf.gz.log"
    threads : 1
    resources:
        mem_mb = 100000
    shell:
        "{params.gatk} --java-options \"-Xmx40g  -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp \" FilterMutectCalls"
        " -V {input.Mutect2_vcf}"
        " -R {input.index}"
        " --contamination-table {input.contamination_table}"
        " -O {output.VCF} 2> {log}"

