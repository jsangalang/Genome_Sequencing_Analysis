# Snakemake:  WGS/WES analysis pipeline for humain/mouse samples

[![Snakemake](https://img.shields.io/badge/snakemake-=5.23.0-brightgreen.svg)](https://snakemake.github.io)

#### Usage On flamingo

- Step 0. clone github workflow project on flamingo
```
$ ssh username@flamingo.intra.igr.fr
$ git clone https://github.com/jinxin-wang/Genome_Sequencing_Analysis.git
```
- Step 1. create conda envirements 
```
$ conda env create -f META_PRISM_conda_env.txt -n META_PRISM
$ conda env create -f Mouse_env.txt -n Mouse
$ conda env create -f pipeline_GATK_2.1.4_conda_env.txt -n pipeline_GATK_2.1.4_V2
```
- Step 2. deploy workflow

Create a directory for your project, then create soft links to your datasets in the directory DNA_samples. Please make sure 
that the name of the files terminated with _0.fastq.gz, _1.fastq.gz or _2.fastq.gz. Moreover, keep in mind that the softlinks 
should be always accessible on the work nodes.

```
$ cd /mnt/beegfs/scratch/username/yourprojectdir
$ mkdir projectname; cd !$
$ ln -s /workflowpath/workflow .
$ cp /workflowpath/run.sh .
$ mkdir -p DNA_samples; cd !$
$ ln -s /datadir/* .
```
- Step 3. configure workflow

1. If there are both tumor samples and normal samples, then you need create a file **variant_call_list_TvN.tsv** in the project directory. The first column is for the tumor samples, and second column is for the normal samples. Each row is a tumor vs normal pair. The seperator is tab. 

2. If there are both tumor samples and panel of normal samples, then you need create a file **variant_call_list_Tp.tsv** in the project directory.

3. If there are tumor, normal samples and panel of normal samples, then you need create a file **variant_call_list_TvNp.tsv** in the project directory.

4. If there are only tumor samples, then the pipeline will run automatically in tumor only mode. 
```
$ cd /mnt/beegfs/scratch/username/yourprojectdir/projectname
$ emacs -nw variant_call_list_TvN.tsv
$ cat variant_call_list.tsv
tumor_sample_A  normal_sample_A
tumor_sample_B  normal_sample_B
$ emacs -nw run.sh
snakemake -c 'sbatch --cpus-per-task={threads} --mem={resources.mem_mb}M -p {params.queue}' --jobs 20 --rerun-incomplete --config samples=mouse seq_type=WES

```
- Step 4. run workflow
```
$ ./run.sh
```

[Examples of Best Practice](https://snakemake.github.io/snakemake-workflow-catalog/)
