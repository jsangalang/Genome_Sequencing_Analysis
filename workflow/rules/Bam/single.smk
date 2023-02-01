## A rule to map single-end DNA sample using BWA
rule bwa_map_single:
    input:
        config["BWA_INDEX"],
        "DNA_samples_clean/{sample}_0.fastq.gz"
    output:
        temp("bam/{sample}.bam")
    log:
        "logs/bam/{sample}.bam.log"
    params:
        queue = "mediumq",    
        bwa = config["APP_BWA"],
        samtools = config["APP_SAMTOOLS"]
    threads: 48
    resources:
        mem_mb = 51200
    shell:
        "{params.bwa} mem -M -R \"@RG\\tID:bwa\\tSM:{wildcards.sample}\\tPL:ILLUMINA\\tLB:truseq\" -t {threads} {input} | {params.samtools} view -@ {threads} -bS - | {params.samtools} sort -@ {threads} - -o {output} 2> {log}"

