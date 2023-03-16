## Annovar on Haplotype caller 
rule annovar:
    input:
        vcf = "haplotype_caller_filtered/{nsample}_germline_variants_filtered.vcf.gz"
    output:
        avinput = "annovar/{nsample}.avinput",
        txt = "annovar/{nsample}.multianno.txt",
        vcf = "annovar/{nsample}.multianno.vcf",
    params:
        queue    = "mediumq",
        annovar  = config["annovar"]["app"],
        ref      = config["annovar"][config["samples"]]["ref"],
        DB       = config["annovar"][config["samples"]]["DB"],
        protocol = config["annovar"][config["samples"]]["protocols"],
        operation= config["annovar"][config["samples"]]["operations"],
    threads : 4
    resources:
        mem_mb = 10240
    log:
        "logs/annovar/{nsample}.log"
    shell :
        "{params.annovar} {input.vcf} {params.DB} "
        " --thread {threads} --maxgenethread 4 "
        " -buildver hg19 "
        " -out annovar/{wildcards.nsample} "
        " -remove "
        " -protocol {params.protocol} "
        " -operation {params.operation} "
        " -nastring . -vcfinput 2> {log} "

