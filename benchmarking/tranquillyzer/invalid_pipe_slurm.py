import os
import typer

def slurm_header(max_time, node, cpus, job_name):
    lines = ["#!/bin/bash",
             f"#SBATCH -t {max_time}",
             f"#SBATCH -p {node}",
             f"#SBATCH --cpus-per-task={cpus}",
             f"#SBATCH --job-name={job_name}",
             "\nsource /home/ayush.semwal/.bashrc",
             "mamba activate tranquilizer",
             "\nmodule load cuda11.8/toolkit/11.8.0",
             "module load bbc2/samtools",
             "module load bbc2/minimap2", "\n"]
    return lines
    
def commands_list(raw_fastq_dir, output_dir, threads, chunk_size, 
                  model, model_type, whitelist):
    
    commands = [f"( time tranquilizer preprocess {raw_fastq_dir} {output_dir} --threads {threads} ) > {output_dir}/preprocess_fasta.log 2>&1",
                f"( time tranquilizer readlengthdist {output_dir} ) > {output_dir}/read_length_distr.log 2>&1",
                f"( time tranquilizer annotate-reads {model} {output_dir} {whitelist} --model-type {model_type} --chunk-size {chunk_size} --threads {threads} ) > {output_dir}/annotate_reads.log 2>&1"]
    
    return commands

app = typer.Typer(rich_markup_mode="rich")

@app.command()
def generate_tranquilizer_pipe(input_base_dir: str = typer.Option(None, help="Input base directory where the fastq files are stored. \
                                                                  The final fastq dir should have the following structure \
                                                                  input_base_dir/bin/sample_size/simulated_data"), 
                               output_base_dir: str = typer.Option(None, help="Output base directory. It will have the following structure\
                                                                   output_base_dir/bin/sample_size/model_type"), 
                               invalid_types: str = typer.Option(None, help="All the invalid types"),
                               threads: int = typer.Option(None, help="Total CPU threads"),
                               chunk_size: int = typer.Option(None, help="Chunk size for read annotation"), 
                               model: str = typer.Option(None, help="Model name for read annotation"), 
                               model_type: str = typer.Option(None, help="Model type for annotation - REG/HYB/CRF"), 
                               max_time: str = typer.Option(None, help="Maximum time allowed for the pipeline to run on the hpc"), 
                               node: str = typer.Option(None, help="Which node to submit the slurm script to"), 
                               cpus: int = typer.Option(None, help="Number of CPU cores to request"), 
                               job_name: str = typer.Option(None, help="What to name the job before submission to HPC")):
    
    invalid_list = invalid_types.split(",")
    for invalid_type in invalid_list:
        raw_fastq_dir = f'{input_base_dir}/{invalid_type}/simulated_data'
        output_dir = f'{output_base_dir}/{invalid_type}/{model_type}'

        os.makedirs(f"{output_dir}", exist_ok=True)

        whitelist = f'{raw_fastq_dir}/barcodes_cleaned.tsv'

        header = slurm_header(max_time, node, cpus, job_name)
        commands = commands_list(raw_fastq_dir, output_dir, threads, chunk_size,
                                 model, model_type, whitelist)
        file_content = "\n".join(header) + "\n\n".join(commands)
        with open(f'{output_dir}/tranquillyzer_pipeline.sh', 'w+') as f:
            f.write(f'{file_content}')
            f.close()

if __name__ == "__main__":
    app()