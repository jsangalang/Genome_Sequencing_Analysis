## A rule to check mapping coverage with mosdepth
rule mosdepth:
    input:
        bam = "bam/{sample}.nodup.recal.bam" if config["remove_duplicates"] == True else "bam/{sample}.recal.bam",
        bai = "bam/{sample}.nodup.recal.bam.bai" if config["remove_duplicates"] == True else "bam/{sample}.recal.bam.bai"
    output:
        "mapping_QC/mosdepth/{sample}.mosdepth.global.dist.txt"
    log:
        "logs/mosdepth/{sample}.log"
    params:
        queue = "mediumq",
        mosdepth = config["mosdepth"]["app"]
    threads : 8
    resources:
        mem_mb = 20480
    shell:
        "{params.mosdepth} -n --fast-mode -t {threads} --by 5000  mapping_QC/mosdepth/{wildcards.sample} {input.bam}"        
