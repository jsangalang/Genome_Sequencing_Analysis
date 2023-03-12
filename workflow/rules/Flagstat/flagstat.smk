## A rule to check mapping metrics with samtools flagstat
rule samtools_flagstat:
    input:
        bam = "bam/{sample}.nodup.recal.bam" if config["REMOVE_DUPLICATES"] == True else "bam/{sample}.recal.bam",
        bai = "bam/{sample}.nodup.recal.bam.bai" if config["REMOVE_DUPLICATES"] == True else "bam/{sample}.recal.bam.bai"
    output:
        "mapping_QC/flagstat/{sample}_flagstat.txt"
    log:
        "logs/mapping_QC/{sample}.flagstat.log"
    params:
        queue = "mediumq",
        samtools = config["SAMTOOLS"]["APP"]
    threads : 16
    resources:
        mem_mb = 51200
    shell:
        "{params.samtools} flagstat -@ {threads} {input.bam} > {output} 2> {log}"
