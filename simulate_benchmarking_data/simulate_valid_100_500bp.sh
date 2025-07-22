#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -p shen
#SBATCH --cpus-per-task=64
#SBATCH --job-name=simulate_benchmarking_data

source /home/ayush.semwal/.bashrc

mamba activate python_3_10_12

cd /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/

# ( time python3 simulate_reads.py 10x3p_sc /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data gencode.v44.transcripts.fa --num-reads 100000 --threads 64  ) > /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/test_simulate.log 2>&1

# ( time python3 simulate_reads.py 10x3p_sc /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc/100_500bp/5_mil gencode.v44.transcripts.fa --cbc-whitelist-domain /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/3M-february-2018.txt --num-reads 2500000 --threads 64  ) > /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc_100_500bp_5_mil_simulate.log 2>&1
# ( time python3 simulate_reads.py 10x3p_sc /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc/100_500bp/25_mil gencode.v44.transcripts.fa --cbc-whitelist-domain /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/3M-february-2018.txt --num-reads 12500000 --threads 64  ) > /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc_100_500bp_25_mil_simulate.log 2>&1
# ( time python3 simulate_reads.py 10x3p_sc /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc/100_500bp/50_mil gencode.v44.transcripts.fa --cbc-whitelist-domain /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/3M-february-2018.txt --num-reads 25000000 --threads 64  ) > /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc_100_500bp_50_mil_simulate.log 2>&1
( time python3 simulate_reads.py 10x3p_sc /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc/100_500bp/75_mil gencode.v44.transcripts.fa --cbc-whitelist-domain /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/3M-february-2018.txt --num-reads 37500000 --threads 64  ) > /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc_100_500bp_75_mil_simulate.log 2>&1
# ( time python3 simulate_reads.py 10x3p_sc /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc/100_500bp/100_mil gencode.v44.transcripts.fa --cbc-whitelist-domain /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/3M-february-2018.txt --num-reads 50000000 --threads 64  ) > /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc_100_500bp_100_mil_simulate.log 2>&1
# ( time python3 simulate_reads.py 10x3p_sc /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc/100_500bp/125_mil gencode.v44.transcripts.fa --num-reads 62500000 --threads 64  ) > /varidata/research/projects/shen/projects/2024_05_02_LRAnnot/simulate_benchmarking_data/10x3p_sc_100_500bp_125_mil_simulate.log 2>&1
