#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -p shen
#SBATCH --cpus-per-task=64
#SBATCH --job-name=simulate_benchmarking_data

source /home/ayush.semwal/.bashrc

mamba activate python_3_10_12

cd /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/


( time python3 simulate_reads.py 10x3p_sc /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc/500_1000bp/10_mil gencode.v44.transcripts.fa --cbc-whitelist-domain /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/3M-february-2018.txt --num-reads 5000000 --min-cdna 500 --max-cdna 999 --threads 64  ) > /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc_500_1000bp_10_mil_simulate.log 2>&1
( time python3 simulate_reads.py 10x3p_sc /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc/1000_1500bp/10_mil gencode.v44.transcripts.fa --cbc-whitelist-domain /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/3M-february-2018.txt --num-reads 5000000 --min-cdna 1000 --max-cdna 1499 --threads 64  ) > /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc_1000_1500bp_10_mil_simulate.log 2>&1
( time python3 simulate_reads.py 10x3p_sc /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc/1500_2000bp/10_mil gencode.v44.transcripts.fa --cbc-whitelist-domain /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/3M-february-2018.txt --num-reads 5000000 --min-cdna 1500 --max-cdna 1999 --threads 64  ) > /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc_1500_2000bp_10_mil_simulate.log 2>&1
( time python3 simulate_reads.py 10x3p_sc /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc/2000_2500bp/10_mil gencode.v44.transcripts.fa --cbc-whitelist-domain /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/3M-february-2018.txt --num-reads 5000000 --min-cdna 2000 --max-cdna 2499 --threads 64  ) > /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc_2000_2500bp_10_mil_simulate.log 2>&1

