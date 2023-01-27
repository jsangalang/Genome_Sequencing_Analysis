## A rule to clean fastq files with fastp PE

rule fastp_PE:
    input:
        fastq_1="DNA_samples/{sample}_1.fastq.gz",
        fastq_2="DNA_samples/{sample}_2.fastq.gz"
    output:
        fastq_clean_1='DNA_samples_clean/{sample}_1.fastq.gz',
        fastq_clean_2='DNA_samples_clean/{sample}_2.fastq.gz',
        html_report='fastp_reports/{sample}_fastp_report.html',
        json_report='fastp_reports/{sample}_fastp_report.json'
    log:
        "logs/fastp/{sample}_fastp.html.log"
    params:
        queue = "longq",
        fastp = config["APP_FASTP"],
        adapters = config["FASTP_ADAPTERS"]
    threads : 16
    resources:
        mem_mb = 51200
    shell:
        '{params.fastp} --thread {threads} --dont_overwrite -i {input.fastq_1} -o {output.fastq_clean_1} -I {input.fastq_2} -O {output.fastq_clean_2} --compression 9 --adapter_fasta {params.adapters} --trim_poly_g --trim_poly_x --length_required 25 --overrepresentation_analysis  --html {output.html_report} --json {output.json_report} 2> {log}'

