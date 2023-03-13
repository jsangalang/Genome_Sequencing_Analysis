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
```
$ cd /mnt/beegfs/scratch/username/yourprojectdir
$ mkdir projectname
$ cd projectname
$ ln -s /workflowpath/workflow .
$ ln -s /workflowpath/conf .
$ cp /workflowpath/run.sh .
$ mkdir -p DNA_samples
$ cd DNA_samples
$ ln -s /datadir/* .
```
- Step 3. configure workflow
If there are only tumor samples, then the pipeline will run automatically in tumor only mode. 
If there are both tumor samples and normal samples, then you need create a file called variant_call_list.tsv in the project directory. 
```
$ emacs -nw run.sh
```
- Step 4. run workflow
```
$ ./run.sh
```

[Examples of Best Practice](https://snakemake.github.io/snakemake-workflow-catalog/)
