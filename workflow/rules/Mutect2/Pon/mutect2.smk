## A rule to call somatic SNPs and indels via local re-assembly of haplotypes, on tumor versus normal tissu using a Panel of Normal, with GATK Mutect2
## Use list of target from target_interval_list.bed if it is present in the working directory

rule Mutect2_pon:
    input:
        index = config["GENOME_FASTA"],
        tumor_bam = "bam/{tsample}.nodup.recal.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{tsample}.recal.bam",
        norm_bam = "bam/{nsample}.nodup.recal.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{nsample}.recal.bam",
        tumor_bai = "bam/{tsample}.nodup.recal.bam.bai" if config["REMOVE_DUPLICATES"]=="True" else "bam/{tsample}.recal.bam.bai",
        norm_bai = "bam/{nsample}.nodup.recal.bam.bai" if config["REMOVE_DUPLICATES"]=="True" else "bam/{nsample}.recal.bam.bai",
        panel_of_normal = "PoN/{panel_of_normal}.vcf",
        interval = config["MUTECT_INTERVAL_DIR"] + "/{interval}.bed",
        GNOMAD_REF = config["GNOMAD_REF"]
    output:
        VCF = "Mutect2_TvNp_tmp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp_ON_{interval}.vcf.gz",
        INDEX = "Mutect2_TvNp_tmp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp_ON_{interval}.vcf.gz.tbi",
        STATS = "Mutect2_TvNp_tmp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp_ON_{interval}.vcf.gz.stats"
    params:
        queue = "mediumq",
        gatk = config["APP_GATK"]
    log:
        "logs/Mutect2_TvNp_tmp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp_ON_{interval}.vcf.log"
    threads : 24
    resources:
        mem_mb = 51200
    shell:
        "read readGroup_{wildcards.tsample} < <(samtools view -H {input.tumor_bam}| grep \'^@RG\' | awk -F\'SM:\' \'{{split($2,a,\" \"); print a[1]}}\' -);"
        "read readGroup_{wildcards.nsample} < <(samtools view -H {input.norm_bam}| grep \'^@RG\' | awk -F\'SM:\' \'{{split($2,a,\" \"); print a[1]}}\' -);"
        "{params.gatk} --java-options \"-Xmx40g  -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp \" Mutect2"
        " --dont-use-soft-clipped-bases true"
        " --native-pair-hmm-threads {threads} "
        " -L {input.interval}"
        " --reference {input.index} "
        " -pon {input.panel_of_normal}"
        " --germline-resource {input.GNOMAD_REF}"
        " -I {input.tumor_bam}"
        " -I {input.norm_bam}"
        " -tumor $readGroup_{wildcards.tsample}"
        " -normal $readGroup_{wildcards.nsample}"
        " -O {output.VCF} 2> {log}"

rule concatenate_mutect2_pon:
    input:
        vcfs = expand("Mutect2_TvNp_tmp/{{tsample}}_Vs_{{nsample}}_PON_{{panel_of_normal}}_TvNp_ON_{mutect_interval}.vcf.gz", mutect_interval=mutect_intervals),
        panel_of_normal = "PoN/{panel_of_normal}.vcf"
    output:
        concatened_vcf = "Mutect2_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp.vcf.gz",
        vcf_liste = "mutect2_TvNp_tmp_list/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp_mutect2_tmp.list"
    params:
        queue = "shortq",
        gatk = config["APP_GATK"]
    threads : 1
    resources:
        mem_mb = 10000
    log:
        "logs/vcftools/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp.vcf.log"
    shell :
        "ls -1a Mutect2_TvNp_tmp/{wildcards.tsample}_Vs_{wildcards.nsample}_PON_{wildcards.panel_of_normal}_TvNp_ON_*gz > mutect2_TvNp_tmp_list/{wildcards.tsample}_Vs_{wildcards.nsample}_PON_{wildcards.panel_of_normal}_TvNp_mutect2_tmp.list &&"
        "{params.gatk}  --java-options \"-Xmx10g  -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp \" MergeVcfs -I {output.vcf_liste} -O {output.concatened_vcf} 2> {log}"

rule concatenate_mutect2_pon_stats:
    input:
        vcfs = expand("Mutect2_TvNp_tmp/{{tsample}}_Vs_{{nsample}}_PON_{{panel_of_normal}}_TvNp_ON_{mutect_interval}.vcf.gz.stats", mutect_interval=mutect_intervals),
        panel_of_normal = "PoN/{panel_of_normal}.vcf"
    output:
        concatened_stats = "Mutect2_TvNp/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp.vcf.gz.stats",
        stat_liste = "mutect2_TvNp_tmp_list/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp_mutect2_tmp_stats.list"
    params:
        queue = "shortq",
        gatk = config["APP_GATK"]
    threads : 1
    resources:
        mem_mb = 10000
    log:
        "logs/vcftools/{tsample}_Vs_{nsample}_PON_{panel_of_normal}_TvNp.vcf.log"
    shell :
        "ls -1a Mutect2_TvNp_tmp/{wildcards.tsample}_Vs_{wildcards.nsample}_PON_{wildcards.panel_of_normal}_TvNp_ON_*stats > mutect2_TvNp_tmp_list/{wildcards.tsample}_Vs_{wildcards.nsample}_PON_{wildcards.panel_of_normal}_TvNp_mutect2_tmp_stats.list &&"
        "{params.gatk}  --java-options \"-Xmx10g  -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp \" MergeMutectStats --stats {output.stat_liste} -O {output.concatened_stats} 2> {log}"

