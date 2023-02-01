## A rule to do Base Quality Score Recalibration (BQSR) - second pass, GATK BaseRecalibrator
rule base_recalibrator_pass2:
    input:
        index = config["GENOME_FASTA"],
        bam = "bam/{sample}.nodup.recal.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{sample}.recal.bam",
        bai = "bam/{sample}.nodup.recal.bam.bai" if config["REMOVE_DUPLICATES"]=="True" else "bam/{sample}.recal.bam.bai",
        GNOMAD_REF = config["GNOMAD_REF"]
    output:
        "BQSR/{sample}_BQSR_pass2.table"
    params:
        queue = "mediumq",
        target_intervals = TARGET_INTERVAL_BQSR,
        gatk = config["APP_GATK"]
    log:
        "logs/BQSR/{sample}_BQSR_pass2.log"
    threads: 16
    resources:
        mem_mb = 51200
    shell:
        "{params.gatk} --java-options \"-Xmx40g -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp -XX:+UseParallelGC -XX:ParallelGCThreads={threads} \" BaseRecalibrator "
        " -R {input.index}"
        " {params.target_intervals}"
        " --known-sites {input.GNOMAD_REF}"
        " -I {input.bam}"
        " -O {output} 2> {log}"

## A rule to analyse covariate BQSR - GATK AnalyzeCovariates
rule analyze_covariates_bqsr:
    input:
        table1 = "BQSR/{sample}_BQSR_pass1.table",
        table2 = "BQSR/{sample}_BQSR_pass2.table"
    output:
        "BQSR/{sample}_BQSR_report.pdf"
    params:
        queue = "mediumq",
        gatk = config["APP_GATK"]
    log:
        "logs/BQSR/{sample}_AnalyzeCovariates.log"
    threads : 16
    resources:
        mem_mb = 51200
    shell:
        "{params.gatk} --java-options \"-Xmx40g -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp -XX:+UseParallelGC -XX:ParallelGCThreads={threads} \" AnalyzeCovariates"
        " --before-report-file {input.table1}"
        " --after-report-file {input.table2}"
        " --plots-report-file {output} 2> {log}"

