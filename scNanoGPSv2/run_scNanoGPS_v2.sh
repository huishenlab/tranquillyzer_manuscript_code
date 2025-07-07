#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=64
#SBATCH -o scNanoGPS.o
#SBATCH -e scNanoGPS.e

# this is a bash script showing the call framework to scNanoGPS

# NUMBER OF READ # --- only parameter that needs chaning & the -o option
nreads=5_mil # e.g. 5_mil 25_mil 50_mil 75_mil 100_mil, 3p_repeat, 500-1000bp etc.

# set required variables
P_DIR="2025_03_12_tranquilizer_benchmarking/scnanogps/scNanoGPS/"
cd $P_DIR
REF_GENOME="2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.fa"
IND_GENOME="2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.mmi"
GENOME_ANNOTATION="2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.gtf"
ncores=64
ANNOVAR="2025_03_12_tranquilizer_benchmarking/2025_03_12_tranquilizer_benchmarking/scnanogps/annovar"
ANNOVAR_DB="2025_03_12_tranquilizer_benchmarking/2025_03_12_tranquilizer_benchmarking/scnanogps/annovar/hg38db/"
ANNOVAR_GV="hg38"
ANNOVAR_PROTOCOL="refGene,cytoBand,gnomad30_genome,avsnp150,dbnsfp42c,cosmic96_coding,cosmic96_noncoding"
ANNOVAR_OP="gx,r,f,f,f,f,f"
ANNOVAR_XREF="tools/annovar/hg38db/omim/gene_xref.txt"

# initialize conda env
conda init bash
source ~/.bashrc  # Ensures the Conda environment is initialized
conda activate scNanoGPS_v2    

echo "Active Conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"
echo "Active python: $(which python)"
echo "Active python version: $(python -V)"
echo "Active conda version: $(conda --version)"

# Note, all this python script does is create a custom bash script with these parameters...instead I run the individual python tools below
# also, option -d does not work
#~ python run_scNanoGPS.py \
#~ -i ${FASTQ} \
#~ -p 3p \
#~ -t ${ncores} \
#~ --gtf=${GENOME_ANNOTATION} \
#~ --ref_genome=${REF_GENOME} \
#~ --idx_genome=${IND_GENOME} \
#~ --whitelist=${FASTQ}cbc_whitelist.txt \
#~ -d scNanoGPS_res_${n_reads}

#_______________________________________________________________________
# run_scNanoGPS.py generates the following bash commands, included here:
P_DIR="2025_03_12_tranquilizer_benchmarking/scnanogps/scNanoGPS"

# specify the reference
REF_GENOME="2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.fa"
IND_GENOME="2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.mmi"
GENOME_ANNOTATION="2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.gtf"

# specify the FASTQ files
FASTQ=2025_03_12_tranquilizer_benchmarking/data/real/scNanoGPS_paper/${nreads}/ # FASTQ for real datasets
#~ FASTQ=2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc/100_500bp/${nreads}/simulated_data/ # FASTQ for different sized (# of reads) datasets (e.g. nreads=5_mil)
#~ FASTQ=2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc/invalid_reads/${nreads}/simulated_data/ # FASTQ for simulated invalid reads (e.g. nreads = 3p_repeat)
#~ FASTQ=/home/ian.beddows/2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc/${nreads}/10_mil/simulated_data/ # FASTQ for length simulated data (e.g. nreads = 500_1000bp)

# specify the whitelist
# note that no whitelist was used when processing real data
WBC=2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc/100_500bp/${nreads}/simulated_data/barcodes_cleaned.tsv # whitelist for different sized (# of reads) datasets (e.g. nreads=5_mil)
WBC=2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc/invalid_reads/${nreads}/simulated_data/cbc_whitelist.txt # whitelist for simulated invalid reads (e.g. nreads = 3p_repeat)
WBC=2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc/${nreads}/10_mil/simulated_data/cbc_whitelist.txt # whitelist for length simulated data (e.g. nreads = 500_1000bp)


PT_SEQ="TTTTTTTTTTTT"

# start & time the individual tools for dataset 'nreads'
echo "read length profiler start"
time python3 $P_DIR/other_utils/read_length_profiler.py -i $FASTQ -d ${nreads} &> logs/run_read_length_profiler.log.txt &
echo "scanner.py start" >&2
time python3 $P_DIR/scanner.py -t $ncores -i $FASTQ --pT $PT_SEQ --min_read_length=100 -d ${nreads} &> logs/${nreads}_run_scanner.log.txt
echo "assigner.py start" >&2
#~ time python3 $P_DIR/assigner.py -t $ncores --whitelist $WBC --min_read_no=1 -d ${nreads} --tmp_dir=${nreads}_tmp &> logs/${nreads}_run_assigner.log.txt # with whitelist
time python3 $P_DIR/assigner.py -t $ncores --min_read_no=1 -d ${nreads} --tmp_dir=${nreads}_no_whitelist_tmp &> logs/${nreads}_run_assigner.log.txt # no whitelist (use with real data)
echo "curator.py start" >&2
time python3 $P_DIR/curator.py -t $ncores --ref_genome $REF_GENOME --idx_genome $IND_GENOME -d ${nreads} --tmp_dir=${nreads} &> logs/${nreads}_run_curator.log.txt
# now have demultiplexed data; time for subsequent steps is not included
echo "reporter_expression.py start" >&2
time python3 $P_DIR/reporter_expression.py -t $ncores --gtf $GENOME_ANNOTATION -d ${nreads} --tmp_dir=${nreads} &> logs/${nreads}_run_reporter_expression.log.txt
echo "reporter_summary.py start" >&2
time python3 $P_DIR/reporter_summary.py --ref_genome $REF_GENOME --gtf $GENOME_ANNOTATION -d ${nreads} --tmp_dir=${nreads} --qualimap_param "--java-mem-size=200G" &> logs/${nreads}_run_reporter_summary.log.txt

