## A rule to remove duplicated readswith picard, will run only if REMOVE_DUPLICATES is set to True in the configuation file
if config["REMOVE_DUPLICATES"]=="True" :
    rule remove_duplicate:
        input:
            "bam/{sample}.bam"
        output:
            bam = temp("bam/{sample}.nodup.bam"),
            metrics = "remove_duplicate_metrics/{sample}.nodup.metrics"
        threads : 16
        resources:
            mem_mb = 81920
        params:
            queue = "mediumq",
            gatk = config["APP_GATK"]
        log:
            "logs/remove_duplicate_metrics/{sample}.nodup.log"
        shell:
            "{params.gatk} --java-options \"-Xmx75g -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp \" MarkDuplicates --INPUT {input} --REMOVE_DUPLICATES=true --OUTPUT {output.bam} -M {output.metrics} 2> {log}"

## A rule to generate bam index with samtools
rule indexbam_before_recal:
    input:
        bam = "bam/{sample}.nodup.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{sample}.bam"
    output:
        bai = temp("bam/{sample}.nodup.bam.bai") if config["REMOVE_DUPLICATES"]=="True" else temp("bam/{sample}.bam.bai")
    params:
        queue = "shortq",
        samtools = config["APP_SAMTOOLS"]
    threads : 12
    resources:
        mem_mb = 30960
    log:
        "logs/bam/{sample}.bam.bai.log"
    shell:
        "{params.samtools} index -@ {threads} {input} 2> {log}"

## A rule to do Base Quality Score Recalibration (BQSR) - first pass, GATK BaseRecalibrator
rule base_recalibrator_pass1:
    input:
        index = config["GENOME_FASTA"],
        bam = "bam/{sample}.nodup.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{sample}.bam",
        bai = "bam/{sample}.nodup.bam.bai" if config["REMOVE_DUPLICATES"]=="True" else "bam/{sample}.bam.bai",
        GNOMAD_REF = config["GNOMAD_REF"]
    output:
        "BQSR/{sample}_BQSR_pass1.table"
    params:
        queue = "mediumq",
        TARGET_INTERVAL_GATK = TARGET_INTERVAL_BQSR,
        gatk = config["APP_GATK"]
    log:
        "logs/BQSR/{sample}_BQSR_pass1.log"
    threads: 16
    resources:
        mem_mb = 51200
    shell:
        "{params.gatk} --java-options \"-Xmx40g -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp \" BaseRecalibrator "
        " -R {input.index}"
        " {params.TARGET_INTERVAL_GATK}"
        " --known-sites {input.GNOMAD_REF}"
        " -I {input.bam}"
        " -O {output} 2> {log}"

## A rule to do Base Quality Score Recalibration (BQSR) - first pass, GATK ApplyBQSR
rule apply_bqsr_pass1:
    input:
        index = config["GENOME_FASTA"],
        bam = "bam/{sample}.nodup.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{sample}.bam",
        bai = "bam/{sample}.nodup.bam.bai" if config["REMOVE_DUPLICATES"]=="True" else "bam/{sample}.bam.bai",
        recal_file = "BQSR/{sample}_BQSR_pass1.table"
    output:
        temp("bam/{sample}.nodup.recal.beforeReformat.bam") if config["REMOVE_DUPLICATES"]=="True" else temp("bam/{sample}.recal.beforeReformat.bam")
    params:
        queue = "mediumq",
        gatk = config["APP_GATK"]
    log:
        "logs/BQSR/{sample}_ApplyBQSR_pass1.log"
    threads: 16
    resources:
        mem_mb = 51200
    shell:
        "{params.gatk} --java-options \"-Xmx40g  -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp -XX:+UseParallelGC -XX:ParallelGCThreads={threads} \" ApplyBQSR "
        " --create-output-bam-index false"
        " -R {input.index}"
        " --bqsr-recal-file {input.recal_file}"
        " -I {input.bam}"
        " -O {output} 2> {log} "
        #" && rm {input.bam} {input.bai}"

## A rule to reformat bam after Base Quality Score Recalibration (BQSR) - with samtools (save a lot of space)
rule reformat_bam:
    input:
        bam = "bam/{sample}.nodup.recal.beforeReformat.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{sample}.recal.beforeReformat.bam",
    output:
        "bam/{sample}.nodup.recal.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{sample}.recal.bam"
    log:
        "logs/bam_reformat/{sample}_reformat.log"
    params:
        queue = "mediumq",
        samtools = config["APP_SAMTOOLS"]
    threads: 16
    resources:
        mem_mb = 51200
    shell:
        "{params.samtools} view -@ {threads} -b -h -o {output} {input.bam}"
        
## A rule to generate bam index with samtools
rule indexbam_after_recal:
    input:
        bam = "bam/{sample}.nodup.recal.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{sample}.recal.bam"
    output:
        bai = "bam/{sample}.nodup.recal.bam.bai" if config["REMOVE_DUPLICATES"]=="True" else "bam/{sample}.recal.bam.bai"
    params:
        queue = "shortq",
        samtools = config["APP_SAMTOOLS"]
    threads : 16
    resources:
        mem_mb = 51200
    log:
        "logs/bam/{sample}.bam.bai.log"
    shell:
        "{params.samtools} index -@ {threads} {input} 2> {log}"