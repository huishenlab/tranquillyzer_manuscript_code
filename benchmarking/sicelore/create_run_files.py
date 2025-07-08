import argparse
import yaml
import os

def read_config(fname):
    """Read YAML config file for files/directories shared across runs.

    Inputs -
        fname - str
    Outputs -
        dict
    """
    with open(fname, 'r') as fh:
        try:
            out = yaml.safe_load(fh)
        except yaml.YAMLError as err:
            print(err)
            exit(1)

    return out

def read_template_file(fname):
    """Read whole file into one string.

    Inputs -
        fname - str
    Returns -
        str
    """
    with open(fname, 'r') as fh:
        out = fh.read()

    return out

def write_file(data, rundir, fname):
    """Write data to output file (creates output directory if it doesn't exist).

    Inputs -
        data   - str
        rundir - str
        fname  - str
    Returns -
        None
    """
    os.makedirs(rundir, exist_ok=True)
    with open(f'{rundir}/{fname}', 'w') as fh:
        fh.write(data)

    return None

def parse_cli():
    """CLI for script.

    Inputs -
        None
    Returns -
        argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(
        prog='create_run_files.py',
    )

    parser.add_argument(
        '-w', '--whitelist',
        default=None,
        help='path to cell barcode whitelist',
    )

    parser.add_argument(
        'config',
        help='YAML for shared files/directories',
    )

    parser.add_argument(
        'rundir',
        help='directory where created files are placed and all files are run from',
    )

    parser.add_argument(
        'fastqs',
        help='directory with data FASTQs',
    )

    parser.add_argument(
        'prefix',
        help='prefix prepended to some output files',
    )

    return parser.parse_args()

def conda_env(rundir):
    """Write conda env file to output directory.

    Inputs -
        rundir - str
    Returns -
        None
    """
    data = read_template_file('templates/template.conda_env.yaml')

    write_file(
        data,
        rundir,
        'conda_env.yaml'
    )

    return None

def nextflow_pipeline(rundir):
    """Write nextflow pipeline file to output directory.

    Inputs -
        rundir - str
    Returns -
        None
    """
    data = read_template_file('templates/template.main.nf')

    write_file(
        data,
        rundir,
        'main.nf'
    )

    return None

def nextflow_config(rundir, args, config):
    """Create config file for nextflow pipeline

    Inputs -
        rundir - str
        args   - argparse.ArgumentParser
        config - dict
    Returns -
        None
    """
    data = read_template_file('templates/template.nextflow.config')

    data = data.replace('WORK'    , f'{rundir}/{config["workdir"]}')
    data = data.replace('FASTQS'  , args.fastqs)
    data = data.replace('PREFIX'  , args.prefix)
    data = data.replace('RUNDIR'  , rundir)
    data = data.replace('OUTPUT'  , config["outdir"])
    data = data.replace('SICELORE', config["sicelore"])
    data = data.replace('JUNCTION', config["juncbed"])
    data = data.replace('REFFA'   , config["minimapfasta"])
    data = data.replace('REFFLAT' , config["refflat"])
    data = data.replace('CBC_LIST', f'"{args.whitelist}"' if args.whitelist is not None else 'null')

    write_file(
        data,
        rundir,
        'nextflow.config'
    )

    return None

def slurm_submit(rundir, jobname, partition):
    """Create slurm submit script for pipeline.

    Inputs -
        rundir    - str
        jobname   - str
        partition - str
    Returns -
        None
    """
    data = read_template_file('templates/template.submit_job.slurm')

    data = data.replace('JOBNAME', jobname)
    data = data.replace('PARTITION', partition)

    write_file(
        data,
        rundir,
        'submit_job.slurm'
    )

    return None

def main():
    args = parse_cli()
    config = read_config(args.config)

    basedir = os.getcwd()

    conda_env(f'{basedir}/{args.rundir}')
    nextflow_pipeline(f'{basedir}/{args.rundir}')
    nextflow_config(f'{basedir}/{args.rundir}', args, config)
    slurm_submit(f'{basedir}/{args.rundir}', args.prefix, config['partition'])

    return None

if __name__ == '__main__':
    main()
