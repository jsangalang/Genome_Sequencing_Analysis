rule annovar_on_mutect_tumor_only:
    input:
        vcf = "Mutect2_T/{tsample}_tumor_only_twicefiltered_T.vcf.gz"
    output:
        avinput = "annovar_mutect2_T/{tsample}.avinput",
        txt = "annovar_mutect2_T/{tsample}.mm9_multianno.txt",
        vcf = "annovar_mutect2_T/{tsample}.mm9_multianno.vcf"
    params:
        queue = "mediumq",
        annovar = config["APP_ANNOVAR"],
        annovar_db = config["ANNOVAR_DB"]
    threads : 4
    resources:
        mem_mb = 20000
    log:
        "logs/annovar_mutect2_T/{tsample}.log"
    shell :
       "{params.annovar} {input.vcf} {params.annovar_db}  --thread 4 --maxgenethread 4 -buildver mm9 -out annovar_mutect2_T/{wildcards.tsample} -remove -protocol refGene,snp128 -operation g,f -nastring . -vcfinput 2> {log}"


