## Annovar on Haplotype caller 
rule annovar:
    input:
        vcf = "haplotype_caller_filtered/{nsample}_germline_variants_filtered.vcf.gz"
    output:
        avinput = "annovar/{nsample}.avinput",
        txt = "annovar/{nsample}.hg19_multianno.txt",
        vcf = "annovar/{nsample}.hg19_multianno.vcf"
    params:
        queue = "mediumq",
        annovar = config["APP_ANNOVAR"],
        annovar_db = config["ANNOVAR"][config["samples"]]["DB"],
    threads : 4
    resources:
        mem_mb = 8000
    log:
        "logs/annovar/{nsample}.log"
    shell :
        "{params.annovar} {input.vcf} {params.annovar_db}  --thread {threads} --maxgenethread 4 -buildver hg19 -out annovar/{wildcards.nsample} -remove -protocol refGene,avsnp150,gnomad211_genome,gnomad211_exome,gme,mcap,revel,regsnpintron,gerp++gt2,clinvar_20200316,intervar_20180118,popfreq_max_20150413,dbnsfp35a,cosmic70,icgc21 -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput 2> {log}"

