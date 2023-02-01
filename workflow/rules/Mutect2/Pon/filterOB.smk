# A rule to filter VCF on orientation bias, for OxoG and FFPE, from mutect tumor Vs normal with PoN 
rule Filter_By_Orientation_Bias_pon:
    input :
        Mutect2_vcf = "Mutect2_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_filtered_TvNp.vcf.gz",
        pre_adapter_detail_metrics = "collect_Sequencing_Artifact_Metrics/{tsample}_artifact.pre_adapter_detail_metrics.txt",
        panel_of_normal = "PoN/{panel_of_normal}.vcf"
    output:
        filtered_vcf = "Mutect2_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_twicefiltered_TvNp.vcf.gz",
        filtered_vcf_index = "Mutect2_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_twicefiltered_TvNp.vcf.gz.tbi"
    params:
        queue = "mediumq",
        gatk = config["APP_GATK"]
    log:
        "logs/filter_Mutect2_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_twicefiltered_TvNp.vcf.gz.log"
    threads : 1
    resources:
        mem_mb = 100000
    shell:
        "{params.gatk} --java-options \"-Xmx40g  -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp \" FilterByOrientationBias"
        " -V {input.Mutect2_vcf}"
        " -AM G/T -AM C/T"
        " -P {input.pre_adapter_detail_metrics}"
        " -O {output.filtered_vcf} 2> {log}"  

