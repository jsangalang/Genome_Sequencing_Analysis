## [The SAM/BAM Format Specification](https://samtools.github.io/hts-specs/SAMtags.pdf)

[Learn bam file format](https://bookdown.org/content/24942ad6-9ed7-44e9-b214-1ea8ba9f0224/learning-the-bam-format.html)

[sam/bam file format explanation](https://genome.sph.umich.edu/wiki/SAM)

[SAM/BAM/CRAM Format](https://learn.gencore.bio.nyu.edu/ngs-file-formats/sambam-format/)

Single-end or paired-end DNA sample mapping using BWA

BWA (Burrows-Wheeler Aligner) is a widely used software package for aligning next-generation sequencing (NGS) data to a reference genome , it supports both single-end and paired-end reads and provides different algorithms tailored for different types of sequencing data
Before mapping the reads, you need to index the reference genome using BWA. This step prepares the reference genome for efficient alignment. 
```
    input:
        fastq = ["DNA_samples_clean/{sample}_1.fastq.gz", "DNA_samples_clean/{sample}_2.fastq.gz"] if config["paired"] == True else "DNA_samples_clean/{sample}_0.fastq.gz",
    output:
        temp("bam/{sample}.bam")
```
The main algorithm implemented in BWA here is :
BWA-MEM: designed for aligning longer reads (e.g., Illumina HiSeq, NovaSeq) to a reference genome.

```
shell:
        "{params.bwa} mem -M -R \"@RG\\tID:bwa\\tSM:{wildcards.sample}\\tPL:ILLUMINA\\tLB:truseq\" -t {threads} {params.index} {input.fastq} | {params.samtools} view -bS - | {params.samtools} sort -@ {threads} - -o {output} 2> {log}"
```
