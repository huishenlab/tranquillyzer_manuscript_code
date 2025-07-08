#!/bin/bash
#SBATCH -J JOBNAME
#SBATCH -o slurm_output.log
#SBATCH --ntasks 1
#SBATCH --time 48:00:00
#SBATCH --mem=10G
#SBATCH --partition=PARTITION

source ~/start_anaconda.sh
conda activate nextflow

cd ${SLURM_SUBMIT_DIR}

nextflow -config nextflow.config run \
    -profile conda \
    main.nf
