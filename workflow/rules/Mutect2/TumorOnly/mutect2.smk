## A rule to call somatic SNPs and indels via local re-assembly of haplotypes, on tumor sample only, with GATK Mutect2
## Use list of target from target_interval_list.bed if it is present in the working directory

## [J. WANG] add " --native-pair-hmm-threads {threads} " to shell
##           raise threads 1 to 24
##           raise mem_mb to 51200

rule Mutect2_tumor_only:
    input:
        index = config["GENOME_FASTA"],
        tumor_bam = "bam/{tsample}.nodup.recal.bam" if config["REMOVE_DUPLICATES"]=="True" else "bam/{tsample}.recal.bam",
        tumor_bai = "bam/{tsample}.nodup.recal.bam.bai" if config["REMOVE_DUPLICATES"]=="True" else "bam/{tsample}.recal.bam.bai",
        interval = config["MUTECT_INTERVAL_DIR"] + "/{interval}.bed",
        GNOMAD_REF = config["GNOMAD_REF"]
    output:
        VCF = "Mutect2_T_tmp/{tsample}_tumor_only_T_ON_{interval}.vcf.gz",
        INDEX = "Mutect2_T_tmp/{tsample}_tumor_only_T_ON_{interval}.vcf.gz.tbi",
        STATS = "Mutect2_T_tmp/{tsample}_tumor_only_T_ON_{interval}.vcf.gz.stats"
    params:
        queue = "mediumq",
        gatk = config["APP_GATK"]
    log:
        "logs/Mutect2_T_tmp/{tsample}_tumor_only_T_ON_{interval}.vcf.log"
    threads : 24
    resources:
        mem_mb = 51200
    shell:
        "read readGroup_{wildcards.tsample} < <(samtools view -H {input.tumor_bam}| grep \'^@RG\' | awk -F\'SM:\' \'{{split($2,a,\" \"); print a[1]}}\' -);"
        "{params.gatk} --java-options \"-Xmx40g  -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp \" Mutect2"
        " --native-pair-hmm-threads {threads} "
        " --dont-use-soft-clipped-bases true"
        " -L {input.interval}"
        " --reference {input.index} "
        " --germline-resource {input.GNOMAD_REF}"
        " -I {input.tumor_bam}"
        " -tumor $readGroup_{wildcards.tsample}"
        " -O {output.VCF} 2> {log}" 
        
rule concatenate_mutect2_tumor_only:
    input:
        vcfs = expand("Mutect2_T_tmp/{{tsample}}_tumor_only_T_ON_{mutect_interval}.vcf.gz", mutect_interval=mutect_intervals)
    output:
        concatened_vcf = "Mutect2_T/{tsample}_tumor_only_T.vcf.gz",
        vcf_liste = "mutect2_T_tmp_list/{tsample}_tumor_only_T_mutect2_tmp.list"
    params:
        queue = "shortq",
        gatk = config["APP_GATK"]
    threads : 1
    resources:
        mem_mb = 10000
    log:
        "logs/vcftools/{tsample}_tumor_only_T.vcf.log"
    shell :
        "ls -1a Mutect2_T_tmp/{wildcards.tsample}_tumor_only_T_ON_*gz > mutect2_T_tmp_list/{wildcards.tsample}_tumor_only_T_mutect2_tmp.list &&"
        "{params.gatk}  --java-options \"-Xmx10g  -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp \" MergeVcfs -I {output.vcf_liste} -O {output.concatened_vcf} 2> {log}"

rule concatenate_mutect2_tumor_only_stats:
    input:
        vcfs = expand("Mutect2_T_tmp/{{tsample}}_tumor_only_T_ON_{mutect_interval}.vcf.gz.stats", mutect_interval=mutect_intervals)
    output:
        concatened_stats = "Mutect2_T/{tsample}_tumor_only_T.vcf.gz.stats",
        stat_liste = "mutect2_T_tmp_list/{tsample}_tumor_only_T_mutect2_tmp_stats.list"
    params:
        queue = "shortq",
        gatk = config["APP_GATK"]
    threads : 1
    resources:
        mem_mb = 10000
    log:
        "logs/vcftools/{tsample}_tumor_only_T.vcf.log"
    shell :
        "ls -1a Mutect2_T_tmp/{wildcards.tsample}_tumor_only_T_ON_*stats > mutect2_T_tmp_list/{wildcards.tsample}_tumor_only_T_mutect2_tmp_stats.list &&"
        "{params.gatk}  --java-options \"-Xmx10g  -Djava.io.tmpdir=/mnt/beegfs/userdata/$USER/tmp \" MergeMutectStats --stats {output.stat_liste} -O {output.concatened_stats} 2> {log}"

## A rule to filter variant call, from mutect tumor only
rule filter_mutect_calls_tumor_only:
    input :
        Mutect2_vcf = "Mutect2_T/{tsample}_tumor_only_T.vcf.gz",
        Mutect2_stats = "Mutect2_T/{tsample}_tumor_only_T.vcf.gz.stats",
        contamination_table = "cross_sample_contamination/{tsample}_calculatecontamination.table",
        MUTECT_FILTER_REF = config["MUTECT_FILTER_REF"],
        index = config["GENOME_FASTA"]
    output:
        VCF = "Mutect2_T/{tsample}_tumor_only_filtered_T.vcf.gz",
        INDEX = "Mutect2_T/{tsample}_tumor_only_filtered_T.vcf.gz.tbi"
    params:
        queue = "mediumq",
        gatk = config["APP_GATK"]
    log:
        "logs/filter_Mutect2_T/{tsample}_tumor_only_filtered_T.vcf.gz.log"
    threads : 1
    resources:
        mem_mb = 100000
    shell:
        "{params.gatk} --java-options \"-Xmx40g  -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp \" FilterMutectCalls"
        " -V {input.Mutect2_vcf}"
        " -R {input.index}"
        " --contamination-table {input.contamination_table}"
        " -O {output.VCF} 2> {log}"
        
## A rule to filter VCF on orientation bias, for OxoG and FFPE, from mutect tumor only 
rule Filter_By_Orientation_Bias_tumor_only:
    input :
        Mutect2_vcf = "Mutect2_T/{tsample}_tumor_only_filtered_T.vcf.gz",
        pre_adapter_detail_metrics = "collect_Sequencing_Artifact_Metrics/{tsample}_artifact.pre_adapter_detail_metrics.txt"
    output:
        filtered_vcf = "Mutect2_T/{tsample}_tumor_only_twicefiltered_T.vcf.gz",
        filtered_vcf_index = "Mutect2_T/{tsample}_tumor_only_twicefiltered_T.vcf.gz.tbi"
    params:
        queue = "mediumq",
        gatk = config["APP_GATK"]
    log:
        "logs/filter_Mutect2_T/{tsample}_tumor_only_twicefiltered_T.vcf.gz.log"
    threads : 1
    resources:
        mem_mb = 100000
    shell:
        "{params.gatk} --java-options \"-Xmx40g  -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp \" FilterByOrientationBias"
        " -V {input.Mutect2_vcf}"
        " -AM G/T -AM C/T"
        " -P {input.pre_adapter_detail_metrics}"
        " -O {output.filtered_vcf} 2> {log}"     
