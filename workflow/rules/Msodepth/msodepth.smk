## A rule to check mapping coverage with mosdepth
rule mosdepth:
    input:
        bam = "bam/{sample}.nodup.recal.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{sample}.recal.bam",
        bai = "bam/{sample}.nodup.recal.bam.bai" if config["REMOVE_DUPLICATES"]=="True" else "bam/{sample}.recal.bam.bai"
    output:
        "mapping_QC/mosdepth/{sample}.mosdepth.global.dist.txt"
    log:
        "logs/mosdepth/{sample}.log"
    params:
        queue = "mediumq",
        mosdepth = config["APP_MOSDEPTH"]
    threads : 4
    resources:
        mem_mb = 10000
    shell:
        "{params.mosdepth} -n --fast-mode -t {threads} --by 5000  mapping_QC/mosdepth/{wildcards.sample} {input.bam}"        
