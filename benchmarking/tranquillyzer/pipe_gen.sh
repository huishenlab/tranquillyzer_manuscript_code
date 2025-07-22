( time python tranquilizer_pipe_slurm.py --input-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc \
--output-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc \
--bins 100_500bp --sample-sizes 5_mil,25_mil,50_mil,75_mil,100_mil --threads 64 --chunk-size 300000 --model 10x3p_sc_ont_011 \
--model-type HYB --ref-fasta /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.fa \
--dedup-threshold 2 --max-time 48:00:00 --node shen --cpus 64 --job-name tranquillyzer_bench ) > /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc/HYB_pipe_gen.log 2>&1

( time python tranquilizer_pipe_slurm.py --input-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc \
--output-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc \
--bins 100_500bp --sample-sizes 5_mil,25_mil,50_mil,75_mil,100_mil --threads 64 --chunk-size 300000 --model 10x3p_sc_ont_011 \
--model-type CRF --ref-fasta /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.fa \
--dedup-threshold 2 --max-time 48:00:00 --node shen --cpus 64 --job-name tranquillyzer_bench ) > /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc/HYB_pipe_gen.log 2>&1

# ( time python tranquilizer_pipe_slurm.py --input-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc \
# --output-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc \
# --bins 1000_1500bp,1500_2000bp --sample-sizes 5_mil,25_mil,50_mil,75_mil,100_mil --threads 64 --chunk-size 300000 --model 10x3p_sc_ont_010 \
# --model-type HYB --ref-fasta /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.fa \
# --dedup-threshold 2 --max-time 48:00:00 --node shen --cpus 64 --job-name tranquillyzer_bench ) > /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc/HYB_pipe_gen.log 2>&1

# ( time python tranquilizer_pipe_slurm.py --input-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc \
# --output-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc \
# --bins 1000_1500bp,1500_2000bp --sample-sizes 5_mil,25_mil,50_mil,75_mil,100_mil --threads 64 --chunk-size 300000 --model 10x3p_sc_ont_010 \
# --model-type CRF --ref-fasta /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.fa \
# --dedup-threshold 2 --max-time 48:00:00 --node shen --cpus 64 --job-name tranquillyzer_bench ) > /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc/HYB_pipe_gen.log 2>&1

# ( time python invalid_pipe_slurm.py --input-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc/invalid_reads \
# --output-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc/invalid_reads \
# --invalid-types 3p_repeat,concat_2x,concat_2x_rev_comp,fwd_fwd_rev,fwd_rev_rev,rev_fwd_rev,5p_repeat,concat_2x_rev,concat_3x,fwd_rev_fwd,rev_fwd_fwd,rev_rev_fwd \
# --threads 64 --chunk-size 300000 --model 10x3p_sc_ont_011 \
# --model-type HYB --max-time 48:00:00 --node shen --cpus 64 --job-name tranquillyzer_bench ) > /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc/HYB_pipe_gen.log 2>&1

# ( time python invalid_pipe_slurm.py --input-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc/invalid_reads \
# --output-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc/invalid_reads \
# --invalid-types 3p_repeat,concat_2x,concat_2x_rev_comp,fwd_fwd_rev,fwd_rev_rev,rev_fwd_rev,5p_repeat,concat_2x_rev,concat_3x,fwd_rev_fwd,rev_fwd_fwd,rev_rev_fwd \
# --threads 64 --chunk-size 300000 --model 10x3p_sc_ont_011 \
# --model-type CRF --max-time 48:00:00 --node shen --cpus 64 --job-name tranquillyzer_bench ) > /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc/HYB_pipe_gen.log 2>&1

# ( time python tranquilizer_pipe_slurm.py --input-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc \
# --output-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc \
# --bins 500_1000bp --sample-sizes 10_mil --threads 64 --chunk-size 300000 --model 10x3p_sc_ont_011 \
# --model-type HYB --ref-fasta /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.fa \
# --dedup-threshold 2 --max-time 48:00:00 --node shen --cpus 64 --job-name tranquillyzer_bench ) > /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc/HYB_pipe_gen.log 2>&1

# ( time python tranquilizer_pipe_slurm.py --input-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc \
# --output-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc \
# --bins 500_1000bp --sample-sizes 10_mil --threads 64 --chunk-size 300000 --model 10x3p_sc_ont_011 \
# --model-type CRF --ref-fasta /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.fa \
# --dedup-threshold 2 --max-time 48:00:00 --node shen --cpus 64 --job-name tranquillyzer_bench ) > /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc/HYB_pipe_gen.log 2>&1             

# ( time python tranquilizer_pipe_slurm.py --input-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc \
# --output-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc \
# --bins 1000_1500bp --sample-sizes 10_mil --threads 64 --chunk-size 300000 --model 10x3p_sc_ont_011 \
# --model-type HYB --ref-fasta /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.fa \
# --dedup-threshold 2 --max-time 48:00:00 --node shen --cpus 64 --job-name tranquillyzer_bench ) > /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc/HYB_pipe_gen.log 2>&1

# ( time python tranquilizer_pipe_slurm.py --input-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc \
# --output-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc \
# --bins 1000_1500bp --sample-sizes 10_mil --threads 64 --chunk-size 300000 --model 10x3p_sc_ont_011 \
# --model-type CRF --ref-fasta /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.fa \
# --dedup-threshold 2 --max-time 48:00:00 --node shen --cpus 64 --job-name tranquillyzer_bench ) > /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc/CRF_pipe_gen.log 2>&1             
                   
# ( time python tranquilizer_pipe_slurm.py --input-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc \
# --output-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc \
# --bins 1500_2000bp --sample-sizes 10_mil --threads 64 --chunk-size 300000 --model 10x3p_sc_ont_011 \
# --model-type HYB --ref-fasta /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.fa \
# --dedup-threshold 2 --max-time 48:00:00 --node shen --cpus 64 --job-name tranquillyzer_bench ) > /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc/HYB_pipe_gen.log 2>&1

# ( time python tranquilizer_pipe_slurm.py --input-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc \
# --output-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc \
# --bins 1500_2000bp --sample-sizes 10_mil --threads 64 --chunk-size 300000 --model 10x3p_sc_ont_011 \
# --model-type CRF --ref-fasta /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.fa \
# --dedup-threshold 2 --max-time 48:00:00 --node shen --cpus 64 --job-name tranquillyzer_bench ) > /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc/CRF_pipe_gen.log 2>&1             

# ( time python tranquilizer_pipe_slurm.py --input-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc \
# --output-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc \
# --bins 2000_2500bp --sample-sizes 10_mil --threads 64 --chunk-size 300000 --model 10x3p_sc_ont_011 \
# --model-type HYB --ref-fasta /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.fa \
# --dedup-threshold 2 --max-time 48:00:00 --node shen --cpus 64 --job-name tranquillyzer_bench ) > /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc/HYB_pipe_gen.log 2>&1

# ( time python tranquilizer_pipe_slurm.py --input-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/data/simulated/10x3p_sc \
# --output-base-dir /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc \
# --bins 2000_2500bp --sample-sizes 10_mil --threads 64 --chunk-size 300000 --model 10x3p_sc_ont_011 \
# --model-type CRF --ref-fasta /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/references/hg38_gencode.fa \
# --dedup-threshold 2 --max-time 48:00:00 --node shen --cpus 64 --job-name tranquillyzer_bench ) > /varidata/research/projects/shen/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc/CRF_pipe_gen.log 2>&1             
                        