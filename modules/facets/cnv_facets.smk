rule cnv_facets:
    input:
        gnomad_ref = config["GNOMAD_REF"],
        tumor_bam = "bam/{tsample}.nodup.recal.bam",
        normal_bam = "bam/{nsample}.nodup.recal.bam",
        tumor_bai = "bam/{tsample}.nodup.recal.bam.bai",
        normal_bai = "bam/{nsample}.nodup.recal.bam.bai"
    output:
        csv = "cnv_facets/{tsample}_Vs_{nsample}.vcf.gz",
        csv_index = "cnv_facets/{tsample}_Vs_{nsample}.vcf.gz.tbi"
    log:
        "logs/facets/{tsample}_Vs_{nsample}_facets.log"
    params:
        queue = "mediumq",
        cnv_facet = config["APP_CNV_FACETS"],
        cval = config["CNV_FACETS_CVAL"],
        ref = config["CNV_FACETS_REF"],
        out_pattern = "{tsample}_Vs_{nsample}"
    threads : 16
    resources:
        mem_mb = 102400
    shell:
        'cd cnv_facets; {params.cnv_facet} --snp-nprocs {threads} --gbuild hg19 --cval {params.cval} --snp-vcf {input.gnomad_ref} -t ../{input.tumor_bam} -n ../{input.normal_bam} --out {params.out_pattern}'
 
