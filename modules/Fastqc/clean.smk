## A rule to generate fastq quality control on cleaned fastq PE

rule fastqc_clean:
    input:
        fastq_clean="DNA_samples_clean/{sample}.fastq.gz",
    output:
        'fastq_QC_clean/clean_{sample}_fastqc.html',
        'fastq_QC_clean/clean_{sample}_fastqc.zip'
    log:
        "logs/fastq_QC_clean/{sample}_fastqc.html.log"
    params:    
        queue = "shortq",
        fastqc = config["APP_FASTQC"],
        adapters = config["FASTQC_ADAPTERS"]
    threads : 16
    resources:
        mem_mb = 51200
    shell:
        '{params.fastqc} -t {threads} -a {params.adapters} -o fastq_QC_clean/ {input} 2> {log}'