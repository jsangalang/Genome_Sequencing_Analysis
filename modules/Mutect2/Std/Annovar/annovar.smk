rule annovar_on_mutect:
    input:
        vcf = "Mutect2_TvN/{tsample}_Vs_{nsample}_twicefiltered_TvN.vcf.gz"
    output:
        avinput = "annovar_mutect2_TvN/{tsample}_Vs_{nsample}.avinput",
        txt = "annovar_mutect2_TvN/{tsample}_Vs_{nsample}.mm9_multianno.txt",
        vcf = "annovar_mutect2_TvN/{tsample}_Vs_{nsample}.mm9_multianno.vcf"
    params:
        queue = "mediumq",
        annovar = config["APP_ANNOVAR"],
        annovar_db = config["ANNOVAR_DB"]
    threads : 4
    resources:
        mem_mb = 20000
    log:
        "logs/annovar/{tsample}_Vs_{nsample}.log"
    shell :
       "{params.annovar} {input.vcf} {params.annovar_db}  --thread 4 --maxgenethread 4 -buildver mm9 -out annovar_mutect2_TvN/{wildcards.tsample}_Vs_{wildcards.nsample} -remove -protocol refGene,snp128 -operation g,f -nastring . -vcfinput 2> {log}"

