## A rule to generate fastq quality control on raw fastq

rule fastqc_raw:
    input:
        #lambda wildcards: [os.path.join(FASTQ_DIR[i], x + '.fastq.gz') for i,x in enumerate(FASTQ_SAMPLES) if x == wildcards.fastq_sample]
        fastq='DNA_samples/{fastq_sample}.fastq.gz'
    output:
        'fastq_QC_raw/{fastq_sample}_fastqc.html',
        'fastq_QC_raw/{fastq_sample}_fastqc.zip'
    log:
        "logs/fastq_QC_raw/{fastq_sample}_fastqc.html.log"
    params:    
        queue = "shortq",
        fastqc = config["APP_FASTQC"],
        adapters = config["FASTQC_ADAPTERS"]
    threads : 16
    resources:
        mem_mb = 51200
    shell:
        '{params.fastqc} -t {threads} -a {params.adapters} -o fastq_QC_raw/ {input.fastq} 2> {log}'
