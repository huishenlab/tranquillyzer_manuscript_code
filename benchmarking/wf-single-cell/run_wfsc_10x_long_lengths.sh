#!/bin/bash
#SBATCH --mail-user=*@vai.org
#SBATCH --mail-type=end,fail
#SBATCH --job-name=wfsc_long_reads
#SBATCH --mem=256G
#SBATCH --output wfsc_run_10x_long_reads.o
#SBATCH --error wfsc_run_10x_long_reads.e
#SBATCH -t 3-00:00:00
#SBATCH --cpus-per-task 64
#SBATCH --array=1-4

echo $(date)

module load bbc2/nextflow/nextflow-23.10.1.5891

DIRS=long_reads.txt
DIR=$(sed -n ${SLURM_ARRAY_TASK_ID}p $DIRS)
FASTQ=$DIR/10_mil/simulated_data
REF_GENOME_DIR=refdata-gex-GRCh38-2024-A
KIT=3prime:v3
THREADS=64
NAME=$(basename $DIR)
OUT_DIR=simulated_data_10x_${NAME}

mkdir -p $OUT_DIR
export SINGULARITY_CACHEDIR=nextflow_singularity_storage
export NXF_SINGULARITY_CACHEDIR=nextflow_singularity_storage
export SINGULARITY_TMPDIR=nextflow_singularity_storage
export NXF_HOME=sim_data_work
NEXTFLOW_CONFIG=sim_data.config
WORKSPACE=sim_data_work
VERSION=v3.1.0-g4580fc3

mkdir -p $WORKSPACE

nextflow run \
        -c $NEXTFLOW_CONFIG \
        -w $WORKSPACE \
	-r $VERSION \
        epi2me-labs/wf-single-cell \
        --fastq $FASTQ \
	--expected_cells 1000 \
	--matrix_min_cells 0 \
	--matrix_min_genes 0 \
        --kit $KIT \
        --ref_genome_dir $REF_GENOME_DIR \
        -profile singularity \
        --out_dir $OUT_DIR \
        --threads $THREADS \
	--sample simulated_data \
	-with-trace

echo $(date)
