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
$ conda env create -f META_PRISM_conda_env.txt
$ conda env create -f Mouse_env.txt
$ conda env create -f pipeline_GATK_2.1.4_conda_env.txt
```
- Step 2. deploy workflow
```
$ cd /mnt/beegfs/scratch/username/yourprojectdir
$ mkdir projectname
$ cd projectname
$ ln -s /workflowpath/workflow .
$ ln -s /workflowpath/conf .
$ cp /workflowpath/run.sh .
```
- Step 3. configure workflow
```
$ emacs -nw run.sh
```
- Step 4. run workflow
```
$ ./run.sh
```

[Examples of Best Practice](https://snakemake.github.io/snakemake-workflow-catalog/)
