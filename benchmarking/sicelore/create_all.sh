# Path to base directory with all FASTQs
# (EDIT WITH YOUR PATH)
FQ_BASE=/path/to/fastq/base/dir

# Which sets of files to create
vary_reads=1
vary_length=1
invalid=1
real=1

# Fixed read length, vary number of reads
if [ $vary_reads -eq 1 ]; then
    python create_run_files.py \
        --whitelist ${FQ_BASE}/simulated/10x3p_sc/100_500bp/5_mil/simulated_data/cbc_whitelist.txt \
        config.yaml \
        100_500bp_005M \
        ${FQ_BASE}/simulated/10x3p_sc/100_500bp/5_mil/simulated_data \
        100_500bp_005M

    python create_run_files.py \
        --whitelist ${FQ_BASE}/simulated/10x3p_sc/100_500bp/25_mil/simulated_data/cbc_whitelist.txt \
        config.yaml \
        100_500bp_025M \
        ${FQ_BASE}/simulated/10x3p_sc/100_500bp/25_mil/simulated_data \
        100_500bp_025M

    python create_run_files.py \
        --whitelist ${FQ_BASE}/simulated/10x3p_sc/100_500bp/50_mil/simulated_data/cbc_whitelist.txt \
        config.yaml \
        100_500bp_050M \
        ${FQ_BASE}/simulated/10x3p_sc/100_500bp/50_mil/simulated_data \
        100_500bp_050M

    python create_run_files.py \
        --whitelist ${FQ_BASE}/simulated/10x3p_sc/100_500bp/75_mil/simulated_data/cbc_whitelist.txt \
        config.yaml \
        100_500bp_075M \
        ${FQ_BASE}/simulated/10x3p_sc/100_500bp/75_mil/simulated_data \
        100_500bp_075M

    python create_run_files.py \
        --whitelist ${FQ_BASE}/simulated/10x3p_sc/100_500bp/100_mil/simulated_data/cbc_whitelist.txt \
        config.yaml \
        100_500bp_100M \
        ${FQ_BASE}/simulated/10x3p_sc/100_500bp/100_mil/simulated_data \
        100_500bp_100M
fi

if [ $vary_length -eq 1 ]; then
    python create_run_files.py \
        --whitelist ${FQ_BASE}/simulated/10x3p_sc/500_1000bp/10_mil/simulated_data/cbc_whitelist.txt \
        config.yaml \
        500_1000bp_10M \
        ${FQ_BASE}/simulated/10x3p_sc/500_1000bp/10_mil/simulated_data \
        500_1000bp_10M

    python create_run_files.py \
        --whitelist ${FQ_BASE}/simulated/10x3p_sc/1000_1500bp/10_mil/simulated_data/cbc_whitelist.txt \
        config.yaml \
        1000_1500bp_10M \
        ${FQ_BASE}/simulated/10x3p_sc/1000_1500bp/10_mil/simulated_data \
        1000_1500bp_10M

    python create_run_files.py \
        --whitelist ${FQ_BASE}/simulated/10x3p_sc/1500_2000bp/10_mil/simulated_data/cbc_whitelist.txt \
        config.yaml \
        1500_2000bp_10M \
        ${FQ_BASE}/simulated/10x3p_sc/1500_2000bp/10_mil/simulated_data \
        1500_2000bp_10M

    python create_run_files.py \
        --whitelist ${FQ_BASE}/simulated/10x3p_sc/2000_2500bp/10_mil/simulated_data/cbc_whitelist.txt \
        config.yaml \
        2000_2500bp_10M \
        ${FQ_BASE}/simulated/10x3p_sc/2000_2500bp/10_mil/simulated_data \
        2000_2500bp_10M
fi

# Invalid reads
if [ $invalid -eq 1 ]; then
    BASE_INVALID=${FQ_BASE}/simulated/10x3p_sc/invalid_reads

    for d in `ls ${BASE_INVALID}`; do
        # d   -> use directory name for invalid structure as prefix and output directory name
        # dir -> location of cbc whitelist and fastqs
        dir=${BASE_INVALID}/${d}/simulated_data

        python create_run_files.py \
            --whitelist ${dir}/cbc_whitelist.txt \
            config.yaml \
            ${d} \
            ${dir} \
            ${d}
    done
fi

# Real data
if [ $real -eq 1 ]; then
    python create_run_files.py \
        --whitelist ${FQ_BASE}/real/scNanoGPS_paper/A375_10mil/barcodes_cleaned.tsv \
        config.yaml \
        A375_10mil \
        ${FQ_BASE}/real/scNanoGPS_paper/A375_10mil \
        A375_10mil
fi
