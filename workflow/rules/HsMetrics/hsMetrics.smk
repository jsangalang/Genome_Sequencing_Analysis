## A rule to get mapping stats, with GATK CollectHsMetrics
rule collect_Hs_metrics:
    input:
        bam = "bam/{sample}.nodup.recal.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{sample}.recal.bam",
        bai = "bam/{sample}.nodup.recal.bam.bai" if config["REMOVE_DUPLICATES"]=="True" else "bam/{sample}.recal.bam.bai",
        interval = config["HSMETRICS_INTERVAL"],
        bait = config["HSMETRICS_BAIT"]
    output:
        "mapping_QC/HsMetrics/{sample}_HsMetrics.tsv"
    params:
        queue = "mediumq",
        gatk = config["APP_GATK"]
    log:
        "logs/mapping_QC/{sample}_HsMetrics.log"
    threads : 16
    resources:
        mem_mb = 51200
    shell:
        "{params.gatk} --java-options \"-Xmx40g  -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp \" CollectHsMetrics"
        " -TI {input.interval}"
        " -BI {input.interval}"
        " -I {input.bam}"
        " -O {output} 2> {log}"