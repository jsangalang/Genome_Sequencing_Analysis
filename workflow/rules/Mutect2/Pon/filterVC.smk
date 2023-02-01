## A rule to filter variant call, from mutect tumor Vs normal with PoN
rule filter_mutect_calls_pon:
    input :
        Mutect2_vcf = "Mutect2_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp.vcf.gz",
        Mutect2_stats = "Mutect2_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp.vcf.gz.stats",
        contamination_table = "cross_sample_contamination/{tsample}_calculatecontamination.table",
        panel_of_normal = "PoN/{panel_of_normal}.vcf",
        MUTECT_FILTER_REF = config["MUTECT_FILTER_REF"],
        index = config["GENOME_FASTA"]
    output:
        VCF = "Mutect2_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_filtered_TvNp.vcf.gz",
        INDEX = "Mutect2_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_filtered_TvNp.vcf.gz.tbi"
    params:
        queue = "mediumq",
        gatk = config["APP_GATK"]
    log:
        "logs/filter_Mutect2_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_filtered_TvNp.vcf.gz.log"
    threads : 1
    resources:
        mem_mb = 100000
    shell:
        "{params.gatk} --java-options \"-Xmx40g  -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp \" FilterMutectCalls"
        " -V {input.Mutect2_vcf}"
        " -R {input.index}"
        " --contamination-table {input.contamination_table}"
        " -O {output.VCF} 2> {log}"

