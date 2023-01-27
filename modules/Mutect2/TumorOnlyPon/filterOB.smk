## A rule to filter VCF on orientation bias, for OxoG and FFPE, from mutect tumor only with PoN 
rule Filter_By_Orientation_Bias_tumor_only_pon:
    input :
        Mutect2_vcf = "Mutect2_Tp/{tsample}_PON_{panel_of_normal}_filtered_Tp.vcf.gz",
        pre_adapter_detail_metrics = "collect_Sequencing_Artifact_Metrics/{tsample}_artifact.pre_adapter_detail_metrics.txt",
    output:
        filtered_vcf = "Mutect2_Tp/{tsample}_PON_{panel_of_normal}_twicefiltered_Tp.vcf.gz",
        filtered_vcf_index = "Mutect2_Tp/{tsample}_PON_{panel_of_normal}_twicefiltered_Tp.vcf.gz.tbi"
    params:
        queue = "mediumq",
        gatk = config["APP_GATK"]        
    log:
        "logs/filter_Mutect2_Tp/{tsample}_PON_{panel_of_normal}_twicefiltered_Tp.vcf.gz.log"
    threads : 1
    resources:
        mem_mb = 100000
    shell:
        "{params.gatk} --java-options \"-Xmx40g  -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp \" FilterByOrientationBias"
        " -V {input.Mutect2_vcf}"
        " -AM G/T -AM C/T"
        " -P {input.pre_adapter_detail_metrics}"
        " -O {output.filtered_vcf} 2> {log}"  
        
