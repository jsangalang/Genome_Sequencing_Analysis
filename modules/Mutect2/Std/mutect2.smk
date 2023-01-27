## A rule to call somatic SNPs and indels via local re-assembly of haplotypes, on tumor versus normal tissu, with GATK Mutect2
## Use list of target from target_interval_list.bed if it is present in the working directory
rule Mutect2:
    input:
        index = config["GENOME_FASTA"],
        tumor_bam = "bam/{tsample}.nodup.recal.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{tsample}.recal.bam",
        norm_bam = "bam/{nsample}.nodup.recal.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{nsample}.recal.bam",
        tumor_bai = "bam/{tsample}.nodup.recal.bam.bai" if config["REMOVE_DUPLICATES"]=="True" else "bam/{tsample}.recal.bam.bai",
        norm_bai = "bam/{nsample}.nodup.recal.bam.bai" if config["REMOVE_DUPLICATES"]=="True" else "bam/{nsample}.recal.bam.bai",
        interval = config["MUTECT_INTERVAL_DIR"] + "/{interval}.bed",
        GNOMAD_REF = config["GNOMAD_REF"]
    output:
        VCF = "Mutect2_TvN_tmp/{tsample}_Vs_{nsample}_TvN_ON_{interval}.vcf.gz",
        INDEX = "Mutect2_TvN_tmp/{tsample}_Vs_{nsample}_TvN_ON_{interval}.vcf.gz.tbi",
        STATS = "Mutect2_TvN_tmp/{tsample}_Vs_{nsample}_TvN_ON_{interval}.vcf.gz.stats"
    params:
        queue = "mediumq",
        tumor_group = "{tsample}",
        norm_group = "{nsample}",
        gatk = config["APP_GATK"]
    log:
        "logs/Mutect2_TvN/{tsample}_Vs_{nsample}_TvN_ON_{interval}.vcf.log"
    threads : 24
    resources:
        mem_mb = 51200
    shell: 
        "read readGroup_{wildcards.tsample} < <(samtools view -@ {threads} -H {input.tumor_bam} | grep \'^@RG\' | awk -F\'SM:\' \'{{split($2,a,\" \"); print a[1]}}\' -);"
        "read readGroup_{wildcards.nsample} < <(samtools view -@ {threads} -H {input.norm_bam} | grep \'^@RG\' | awk -F\'SM:\' \'{{split($2,a,\" \"); print a[1]}}\' -);"
        "{params.gatk} --java-options \"-Xmx40g  -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp \" Mutect2"
        " --dont-use-soft-clipped-bases true "
        " --native-pair-hmm-threads {threads} "
        " -L {input.interval}"
        " --reference {input.index} "
        " --germline-resource {input.GNOMAD_REF}"
        " -I {input.tumor_bam}"
        " -I {input.norm_bam}"
        " -tumor $readGroup_{wildcards.tsample}"
        " -normal $readGroup_{wildcards.nsample}"
        " -O {output.VCF} 2> {log}"

## Concatenate mutect2 results
rule concatenate_mutect2:
    input:
        vcfs = expand("Mutect2_TvN_tmp/{{tsample}}_Vs_{{nsample}}_TvN_ON_{mutect_interval}.vcf.gz", mutect_interval=mutect_intervals)
    output:
        concatened_vcf = "Mutect2_TvN/{tsample}_Vs_{nsample}_TvN.vcf.gz",
        vcf_liste = "Mutect2_TvN_tmp_list/{tsample}_Vs_{nsample}_TvN_mutect2_tmp.list"
    params:
        queue = "shortq",
        gatk = config["APP_GATK"]
    threads : 1
    resources:
        mem_mb = 10000
    log:
        "logs/vcftools/{tsample}_Vs_{nsample}_TvN.vcf.log"
    shell :
        "ls -1a Mutect2_TvN_tmp/{wildcards.tsample}_Vs_{wildcards.nsample}_TvN_ON_*gz > Mutect2_TvN_tmp_list/{wildcards.tsample}_Vs_{wildcards.nsample}_TvN_mutect2_tmp.list && "
        "{params.gatk}  --java-options \"-Xmx10g  -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp \" MergeVcfs -I {output.vcf_liste} -O {output.concatened_vcf} 2> {log}"
 
## Concatenation of Mutect2 stats
rule concatenate_mutect2_stats:
    input:
        vcfs = expand("Mutect2_TvN_tmp/{{tsample}}_Vs_{{nsample}}_TvN_ON_{mutect_interval}.vcf.gz.stats", mutect_interval=mutect_intervals)
    output:
        concatened_stats = "Mutect2_TvN/{tsample}_Vs_{nsample}_TvN.vcf.gz.stats",
        stat_liste = "Mutect2_TvN_tmp_list/{tsample}_Vs_{nsample}_TvN_mutect2_tmp_stats.list"
    params:
        queue = "shortq",
        gatk = config["APP_GATK"]
    threads : 1
    resources:
        mem_mb = 10000
    log:
        "logs/vcftools/{tsample}_Vs_{nsample}_TvN.vcf.log"
    shell :
        "ls -1a Mutect2_TvN_tmp/{wildcards.tsample}_Vs_{wildcards.nsample}_TvN_ON_*stats > Mutect2_TvN_tmp_list/{wildcards.tsample}_Vs_{wildcards.nsample}_TvN_mutect2_tmp_stats.list &&"
        "{params.gatk}  --java-options \"-Xmx10g  -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp \" MergeMutectStats --stats {output.stat_liste} -O {output.concatened_stats} 2> {log}"

