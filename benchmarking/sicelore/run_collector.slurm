#!/bin/bash
#SBATCH -J collector
#SBATCH -o run_collector.log
#SBATCH --ntasks 1
#SBATCH --time 48:00:00
#SBATCH --mem=100G
#SBATCH --partition=YOUR_PARTITION_NAME

source ~/start_anaconda.sh

cd ${SLURM_SUBMIT_DIR}

python collect_benchmarking_metrics.py | tee run_collector_output.log
