rule annovar_on_mutect_tumor_only:
    input:
        vcf = "Mutect2_T/{tsample}_tumor_only_twicefiltered_T.vcf.gz"
    output:
        avinput = "annovar_mutect2_T/{tsample}.avinput",
        txt = "annovar_mutect2_T/{tsample}.mm9_multianno.txt",
        vcf = "annovar_mutect2_T/{tsample}.mm9_multianno.vcf"
    params:
        queue = "mediumq",
        annovar = config["annovar"]["app"],
        annovar_db  = config["annovar"][config["samples"]]["DB"],
        annovar_ref = config["annovar"][config["samples"]]["ref"],
    threads : 8
    resources:
        mem_mb = 20480
    log:
        "logs/annovar_mutect2_T/{tsample}.log"
    shell :
       "{params.annovar} {input.vcf} {params.annovar_db} --thread {threads} --maxgenethread 4 -buildver {params.annovar_ref} -out annovar_mutect2_T/{wildcards.tsample} -remove -protocol refGene,snp128 -operation g,f -nastring . -vcfinput 2> {log}"


