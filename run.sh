#!/bin/bash

/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake -c 'sbatch --cpus-per-task={threads} --mem={resources.mem_mb}M -p {params.queue}' --jobs 20 --latency-wait 50 --rerun-incomplete --config samples=human seq_type=WGS

