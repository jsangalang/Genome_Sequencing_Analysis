rule facets_snp_pilleup:
    input:
        GNOMAD_REF = config["GNOMAD_REF"],
        TUMOR_BAM = "bam/{tsample}.nodup.recal.bam",
        NORMAL_BAM = "bam/{nsample}.nodup.recal.bam"
    output:
        CSV = "facets/{tsample}_Vs_{nsample}_facets.csv.gz"
    log:
        "logs/facets/{tsample}_Vs_{nsample}_facets.log"
    params:
        queue = "mediumq",
        snp_pileup = config["APP_FACET_SNP_PILEUP"]
    threads : 1
    resources:
        mem_mb = 10000
    shell:
        '{params.snp_pileup} -g --min-map-quality=55 --min-base-quality=20 --max-depth=200 --min-read-counts=15,15 {input.GNOMAD_REF} {output.CSV} {input.NORMAL_BAM} {input.TUMOR_BAM}'
		
#A rule to draw facets graphs
rule facet_graph:
    input:
        CSV="facets/{tsample}_Vs_{nsample}_facets.csv.gz"
    output:
        PDF =  "facets/{tsample}_Vs_{nsample}_facets_cval500.pdf",
        RDATA = "facets/{tsample}_Vs_{nsample}_facets_cval500.RData"
    log:
        "logs/facets/{tsample}_Vs_{nsample}_facets_graph.log"
    params:
        queue = "mediumq",
        facet_graph = config["APP_FACET_GRAPH"],
        Rscript = config["APP_RSCRIPT"]
    threads : 1
    resources:
        mem_mb = 100000
    shell:
        '{params.Rscript} {params.facet_graph} {input.CSV}'