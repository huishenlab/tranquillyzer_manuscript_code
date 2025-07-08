# Sicelore Benchmarking Code

Scripts used to run Sicelore benchmarking for the tranquillyzer manuscript

## `templates/template.conda_env.yaml`

- Conda environment file used when running pipeline

## `templates/template.main.nf`

- Nextflow pipeline

## `templates/template.nextflow.config`

- Config file used when running pipeline

## `templates/template.submit_job.slurm`

- Submit nextflow pipeline to a SLURM workload manager
- Nextflow version used: 24.10.5

## `collect_benchmarking_metrics.py`

- Collect benchmarking metrics from outputs of pipeline
- Python dependencies
  - `python`: 3.10.16
  - `pandas`: 2.2.2
  - `pysam`: 0.20.0
  - `re`: 2.2.1

## `config.yaml`

- Config file for files/variables shared across pipeline runs

## `create_all.sh`

- Create files for all benchmarking runs

## `create_ref_files.sh`

- Create additional reference files needed for running sicelore
- Dependency versions used:
  - `minimap2`: 2.24
  - `gtfToGenePred`: v479

## `create_run_files.py`

- Create pipeline files using `config.yaml` and template files
- Python dependency versions used
  - `python`: 3.10.16
  - `argparse`: 1.1
  - `yaml`: 6.0.2

## `run_collector.slurm`

- Run benchmarking collect from a SLURM workload manager
