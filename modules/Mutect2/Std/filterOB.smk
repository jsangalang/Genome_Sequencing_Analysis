## A rule to filter VCF on orientation bias, for OxoG and FFPE, from mutect tumor Vs normal 
rule Filter_By_Orientation_Bias:
    input :
        Mutect2_vcf = "Mutect2_TvN/{tsample}_Vs_{nsample}_filtered_TvN.vcf.gz",
        pre_adapter_detail_metrics = "collect_Sequencing_Artifact_Metrics/{tsample}_artifact.pre_adapter_detail_metrics.txt"
    output:
        filtered_vcf = "Mutect2_TvN/{tsample}_Vs_{nsample}_twicefiltered_TvN.vcf.gz",
        filtered_vcf_index = "Mutect2_TvN/{tsample}_Vs_{nsample}_twicefiltered_TvN.vcf.gz.tbi"
    params:
        queue = "mediumq",
        gatk = config["APP_GATK"]
    log:
        "logs/filter_Mutect2_TvN/{tsample}_Vs_{nsample}_twicefiltered_TvN.vcf.gz.log"
    threads : 1
    resources:
        mem_mb = 100000
    shell:
        "{params.gatk} --java-options \"-Xmx40g  -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp \" FilterByOrientationBias"
        " -V {input.Mutect2_vcf}"
        " -AM G/T -AM C/T"
        " -P {input.pre_adapter_detail_metrics}"
        " -O {output.filtered_vcf} 2> {log}"
