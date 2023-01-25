## A rule to clean fastq files with fastp SE

rule fastp_SE:
    input:
        fastq_0="DNA_samples/{sample}_0.fastq.gz",
    output:
        fastq_clean='DNA_samples_clean/{sample}_0.fastq.gz',
        html_report='fastp_reports/{sample}_fastp_report.html',
        json_report='fastp_reports/{sample}_fastp_report.json'
    log:
        "logs/fastp/{sample}_fastp.html.log"
    params:
        queue = "mediumq",
        fastp = config["APP_FASTP"],
        adapters = config["FASTP_ADAPTERS"]
    threads : 16
    resources:
        mem_mb = 51200
    shell:
        '{params.fastp} --thread {threads} --dont_overwrite -i {input.fastq_0} -o {output.fastq_clean} --compression 9 --adapter_fasta {params.adapters} --trim_poly_g --trim_poly_x --length_required 25 --overrepresentation_analysis  --html {output.html_report} --json {output.json_report} 2> {log}'
