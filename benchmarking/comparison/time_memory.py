import re
import os
import polars as pl

def extract_timing_and_memory(logfile):
    with open(logfile, 'r')  as f:
        lines = f.readlines()[-20:]  # Read last 10 lines (safety buffer)
    # result = ()
    for line in lines:
        if "Peak memory usage" in line:
            match = re.search(r"Peak memory usage.*?:\s+([\d.]+)\s+MB", line)
            if match:
                result = (float(match.group(1)),)
        elif line.startswith("real"):
            match = re.search(r"real\s+(\d+)m([\d.]+)s", line)
            if match:
                result = (float(int(match.group(1)) * 60 + float(match.group(2)))/60,) + result
    return result

########### 100_500bp ###########

base_dir = "/home/ayush.semwal/shen-secondary/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc/100_500bp"
sample_sizes = ["5_mil", "25_mil", "50_mil", "75_mil", "100_mil"]
steps = ["preprocess_fasta", "annotate_reads", "align_reads", "dedup"]

approaches = ["HYB", "CRF"]

results = []

for approach in approaches:
    for sample_size in sample_sizes:
        result = [sample_size, approach]
        for step in steps:
            step_res = extract_timing_and_memory(os.path.join(base_dir, sample_size, approach, step +".log"))
            result.append(step_res[0])
            result.append(step_res[1])
        results.append(tuple(result))

summary_df = pl.DataFrame(
    data=results,
    schema=["sample_size", "approach", "preprocess_fasta_time", "preprocess_fasta_mem", "annotate_demux_time", 
    "annotate_demux_mem", "align_time", "align_mem", "dedup_time", "dedup_mem"],
    orient='row'
)

summary_df.write_csv("/home/ayush.semwal/shen-secondary/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/time_and_memory_summary_100_500bp.tsv", separator="\t")

########### longer reads ###########

base_dir = "/home/ayush.semwal/shen-secondary/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/simulated/10x3p_sc"
bin_sizes = ["500_1000bp", "1000_1500bp", "1500_2000bp", "2000_2500bp"]
sample_sizes = ["10_mil"]
steps = ["preprocess_fasta", "annotate_reads", "align_reads", "dedup"]

approaches = ["HYB", "CRF"]

results = []

for approach in approaches:
    for bin_size in  bin_sizes:
        for sample_size in sample_sizes:
            result = [bin_size, sample_size, approach]
            for step in steps:
                step_res = extract_timing_and_memory(os.path.join(base_dir, bin_size, sample_size, approach, step +".log"))
                result.append(step_res[0])
                result.append(step_res[1])
            results.append(tuple(result))

summary_df = pl.DataFrame(
    data=results,
    schema=["bin_size", "sample_size", "approach", "preprocess_fasta_time", "preprocess_fasta_mem", "annotate_demux_time", 
    "annotate_demux_mem", "align_time", "align_mem", "dedup_time", "dedup_mem"],
    orient='row'
)

summary_df.write_csv("/home/ayush.semwal/shen-secondary/projects/2025_03_12_tranquilizer_benchmarking/tranquilizer/time_and_memory_summary_500_2500bp.tsv", separator="\t")