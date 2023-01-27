## A rule to estimate cross-sample contamination using GetPileupSummaries and CalculateContamination, step one GetPileupSummaries
rule get_pileup_summaries:
    input:
        tumor_bam = "bam/{tsample}.nodup.recal.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{tsample}.recal.bam",
        tumor_bai = "bam/{tsample}.nodup.recal.bam.bai" if config["REMOVE_DUPLICATES"]=="True" else "bam/{tsample}.recal.bam.bai",
        MUTECT_FILTER_REF = config["MUTECT_FILTER_REF"]
    output:
        "cross_sample_contamination/{tsample}_getpileupsummaries.table"
    params:
        queue = "mediumq",
        gatk = config["APP_GATK"]
    log:
        "logs/cross_sample_contamination/{tsample}_getpileupsummaries.table.log"
    threads : 1
    resources:
        mem_mb = 15000
    shell:
        "{params.gatk} --java-options \"-Xmx20g  -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp \" GetPileupSummaries"
        " -I {input.tumor_bam}"
        " -L {input.MUTECT_FILTER_REF}"
        " -V {input.MUTECT_FILTER_REF}"
        " -O {output} 2> {log}"
        
##  A rule to estimate cross-sample contamination using GetPileupSummaries and CalculateContamination, step two CalculateContamination
rule calculate_contamination:
    input:
        table = "cross_sample_contamination/{tsample}_getpileupsummaries.table",
    output:
        "cross_sample_contamination/{tsample}_calculatecontamination.table"
    params:
        queue = "mediumq",
        gatk = config["APP_GATK"]
    log:
        "logs/cross_sample_contamination/{tsample}_calculatecontamination.table.log"
    threads : 1
    resources:
        mem_mb = 5000
    shell:
        "{params.gatk} --java-options \"-Xmx40g  -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp \" CalculateContamination"
        " -I {input.table}"
        " -O {output} 2> {log}"

