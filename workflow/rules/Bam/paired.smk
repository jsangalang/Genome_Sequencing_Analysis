## A rule to map paired-end DNA sample using BWA
rule bwa_map_paired:
    input:
        config["BWA_INDEX"],
        "DNA_samples_clean/{sample}_1.fastq.gz",
        "DNA_samples_clean/{sample}_2.fastq.gz"
    output:
        temp("bam/{sample}.bam")
    log:
        "logs/bam/{sample}.bam.log"
    params:
        queue = "mediumq",
        bwa = config["APP_BWA"],
        samtools = config["APP_SAMTOOLS"]
    threads: 16
    resources:
        mem_mb = 51200
    shell:
        "{params.bwa} mem -M -R \"@RG\\tID:bwa\\tSM:{wildcards.sample}\\tPL:ILLUMINA\\tLB:truseq\" -t 16 {input} | {params.samtools} view -@ 16 -bS - | {params.samtools} sort -@ 16 - -o {output} 2> {log}"
