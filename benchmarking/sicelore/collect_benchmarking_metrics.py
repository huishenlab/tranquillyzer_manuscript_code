import subprocess
import pandas as pd
import pysam
import time
import glob
import gzip
import re
import os

OUTPUT_DIR = 'benchmark_results'

def write_readme():
    """Add README describing benchmarking results files

    Inputs -
        None
    Returns -
        None
    """
    with open(f'{OUTPUT_DIR}/README.md', 'w') as fh:
        fh.write('# Barcode Results Files\n')
        fh.write('\n')
        fh.write('- **File name format:** barcodes.<sample_name>.tsv\n')
        fh.write('- **File description:** corrected cell barcodes and UMIs from BAM file generated during UMI assignment step (STEP 3)\n')
        fh.write('- **Column descriptions:**\n')
        fh.write('  - `read_name`: read name extracted from BAM (may not exactly match read name in FASTQ due to parsing sicelore read\n')
        fh.write('  name additions)\n')
        fh.write('  - `cell_barcode`: corrected cell barcode from the BC cell barcode tag in BAM (determined in STEP 1)\n')
        fh.write('  - `umi`: corrected UMI from U8 tag in BAM (determined in STEP 3)\n')
        fh.write('\n')
        fh.write('# Benchmarking Results\n')
        fh.write('\n')
        fh.write('- **File name format:** benchmarking_results.tsv\n')
        fh.write('- **File description:** Statistics related to benchmarking for each sample run\n')
        fh.write('- **Column descriptions:**\n')
        fh.write('  - `dataset`: name of processed sample (generally matches directory name in `data/`)\n')
        fh.write('  - `t_step_1`: time (in seconds) to scan reads for chimeric cDNA, poly(A/T) and adapter sequences, and cell barcodes;\n')
        fh.write('  will assign scanned cell barcodes to a specific cell barcode at this stage - implicitly determines correct read\n')
        fh.write('  architecture at this point\n')
        fh.write('  - `t_step_2`: time (in seconds) to map reads passing step 1 using minimap2\n')
        fh.write('  - `t_step_3`: time (in seconds) to assign UMIs and pull information from read name into SAM auxiliary tags\n')
        fh.write('  - `max_memory_kb`: maximum memory (in kB) used during the pipeline\n')
        fh.write('  - `n_arch_passed`: number of reads in the passed FASTQ from step 1 (i.e., passed architecture and other initial\n')
        fh.write('  filters)\n')
        fh.write('  - `n_arch_failed`: number of reads in the failed FASTQ from step 1 (i.e., did not pass the architecture and other\n')
        fh.write('  initial filters)\n')
        fh.write('  - `n_assigned_bc`: of the `n_arch_passed` reads, the number of reads with a BC extracted from the read name\n')
        fh.write('  - `n_primary`: number of primary alignment reads in BAM from step 2 - the number of alignments from minimap2\n')
        fh.write('  - `n_assigned_umi`: of the `n_primary` reads, the number of reads with a U8 SAM auxiliary tag (essentially the\n')
        fh.write('  number of reads that were primary alignments and have both a BC and a UMI\n')
        fh.write('\n')
        fh.write('# Notes\n')
        fh.write('\n')
        fh.write('- `n_arch_passed + n_arch_failed` may be greater than the number of reads in the original FASTQ due to sicelore\n')
        fh.write('splitting original reads due to chimeras\n')

    return None

def get_latest(t1, t2):
    """Compare two times and return the latest one

    Inputs -
        t1 - tuple (str, int)
        t2 - tuple (str, int)
    Returns -
        tuple (str, int)
    """
    # Separate out filename parts
    pieces_1 = t1[0].split('_')
    pieces_2 = t2[0].split('_')

    # Start by comparing dates
    date_1 = pieces_1[2].split('-')
    date_2 = pieces_2[2].split('-')

    for one, two in zip(date_1, date_2):
        if one != two:
            return t1 if one > two else t2

    # Only compare times if the dates are equal
    time_1 = pieces_1[3].split('-')
    time_2 = pieces_2[3].split('-')

    for one, two in zip(time_1, time_2):
        if one != two:
            return t1 if one > two else t2

    # If somehow we got to here, t1 and t2 are the same, so default to t1
    return t1

def get_pipeline_info(base_name):
    """Extract pipeline run info into a DataFrame

    Inputs -
        base_name - str
    Returns -
        DataFrame
    """
    dir_name = f'{base_name}/final_files/pipeline_info'

    # Retrieve data file for latest running execution_trace* file
    fnames = glob.glob(f'{dir_name}/execution_trace_*')

    fname = None
    if len(fnames) == 1:
        fname = fnames[0]
    else:
        latest = (fnames[0], 0)
        for i in ranges(1, len(fnames)):
            tmp = get_latest(latest, (fnames[i], i))
            if tmp[1] != latest[1]:
                latest = tmp

        fname = latest[0]

    df = pd.read_csv(fname, sep='\t')

    return df

def count_failed_fastq_reads(fname):
    """Count number of lines in a gzipped failed fastq file

    Inputs -
        fname - str
    Returns -
        int
    """
    with gzip.open(fname, 'rb') as fh:
        for i, l in enumerate(fh):
            pass

    try:
        count = i+1
    except UnboundLocalError:
        count = 0

    return count

def read_passed_fastq_file(fname, oname):
    """Count number of passed reads and barcode extracted reads. Output read positions of extracted info

    Inputs -
        fname - str
        oname - str
    Returns -
        tuple, (int, int)
    """
    n_passed = 0
    n_has_bc = 0
    with gzip.open(fname, 'rt') as fh:
        with gzip.open(oname, 'wt', compresslevel=6) as out:
            out.write('name\tstrand\tpolyA_start\tpolyA_end\tadapter_end\ttso_end\tbc_start\tbc_end\n')
            for idx, l in enumerate(fh):
                if idx % 4 == 0:
                    n_passed += 1
                    if n_passed % 100000 == 0:
                        print(f'\t\tprocessed {n_passed} reads')

                    line = l.strip().lstrip('@').split(' ')[0].split('_')

                    name = []
                    strand = None
                    polyA_start = None
                    polyA_end = None
                    adapter_end = None
                    tso_end = None
                    bc_start = None
                    bc_end = None

                    collect_name = True
                    for piece in line:
                        if piece in ['REV', 'FWD']:
                            collect_name = False
                            strand = piece

                        if collect_name:
                            name.append(piece)
                        else:
                            if piece.startswith('PS='):
                                polyA_start = piece.split('=')[1]
                            elif piece.startswith('PE='):
                                polyA_end = piece.split('=')[1]
                            elif piece.startswith('AE='):
                                adapter_end = piece.split('=')[1]
                            elif piece.startswith('T='):
                                tso_end = piece.split('=')[1]
                            elif piece.startswith('bcStart='):
                                bc_start = piece.split('=')[1]
                            elif piece.startswith('bcEnd='):
                                bc_end = piece.split('=')[1]
                            elif piece.startswith('bc='):
                                n_has_bc += 1

                    strand = strand if not None else '.'
                    polyA_start = polyA_start if not None else '.'
                    polyA_end = polyA_end if not None else '.'
                    adapter_end = adapter_end if not None else '.'
                    tso_end = tso_end if not None else '.'
                    bc_start = bc_start if not None else '.'
                    bc_end = bc_end if not None else '.'

                    out.write(f'{"_".join(name)}\t{strand}\t{polyA_start}\t{polyA_end}\t{adapter_end}\t{tso_end}\t{bc_start}\t{bc_end}\n')

    return n_passed, n_has_bc

def get_passed_reads(dir, df):
    """Count number of reads passing/not passing filters

    Inputs -
        dir - str
        df  - DataFrame
    Returns -
        tuple (int, int, int)
    """
    passed_file = f'{dir}/final_files/01.readscan/fastq_pass.fastq.gz'
    n_passed, n_has_bc = read_passed_fastq_file(passed_file, f'{OUTPUT_DIR}/positions.{dir}.tsv.gz')

    hash = df['hash'].iloc[0]
    dir_path = glob.glob(f'{dir}/intermediate_files/{hash}*/passed')
    
    if len(dir_path) > 1:
        print(f'Too many directories to choose from! {dir_path}')

    dir_path = dir_path[0]
    
    n_failed = 0
    failed_files = glob.glob(f'{dir_path}/failed/*.fastq')
    for idx, fname in enumerate(failed_files, start=1):
        if idx % 100 == 0:
            print(f'\t\tprocessed {idx} failed files')
        n_failed += count_failed_fastq_reads(fname)

    return n_passed, n_has_bc, int(n_failed / 4)

def convert_time(time):
    """Convert a time string (Xh Ym Zs) into seconds

    Inputs -
        time - str
    Returns -
        float
    """
    pieces = time.split(' ')

    sec = 0
    for p in pieces:
        m = re.match(r'([0-9.]+)([a-z]+)', p, re.I)
        if not m:
            print(f'Could not parse time: {p} ({time})')
            exit(1)

        number = float(m.group(1))
        indicator = m.group(2)

        if indicator == 'h':
            sec += 3600 * number
        elif indicator == 'm':
            sec += 60 * number
        elif indicator == 's':
            sec += number
        elif indicator == 'ms':
            sec += number / 1000
        else:
            print(f'Unknown time indicator {indicator} ({time})')
            exit(1)

    return sec

def get_timing(df):
    """Get step runtimes (assumes specific order of steps in nextflow trace file)

    Inputs -
        df - DataFrame
    Returns -
        tuple (int, int, int)
    """
    times = list(df['realtime'])

    t1 = convert_time(times[0]) + convert_time(times[1])
    t2 = convert_time(times[2])
    t3 = convert_time(times[3])

    return (t1, t2, t3)

def get_peak_memory(df):
    """Convert memory to consistent units and return peak memory usage

    Inputs -
        df - DataFrame
    Returns -
        float
    """
    mems = list(df['peak_rss'])

    mems_kb = []
    for mem in mems:
        pieces = mem.split(' ')
        if pieces[1] == 'kB':
            mems_kb.append(float(pieces[0]))
        elif pieces[1] == 'MB':
            mems_kb.append(1000 * float(pieces[0]))
        elif pieces[1] == 'GB':
            mems_kb.append(1000000 * float(pieces[0]))
        else:
            print(f'Unknown memory unit found {pieces[1]} ({mem})')
            exit(1)

    return max(mems_kb)

def get_n_primary_alignments(dir):
    """Use samtools to count number of primary alignments

    Inputs -
        dir - str
    Returns -
        int
    """
    bam = f'{dir}/final_files/02.mapping/passed.bam'

    output = subprocess.run(['samtools', 'view', '-c', '-F', '2308', bam], capture_output=True)

    return int(output.stdout)

def get_original_name(name):
    """Get original read name from sicelore mangled name

    Inputs -
        name - str
    Returns -
        str
    """
    pieces = name.split('_')

    out = []
    for piece in pieces:
        start = piece[:2]
        if start in ['RE', 'FW', 'FA', 'FT', 'RA', 'RT']:
            break

        out.append(piece)

    return '_'.join(out)

def process_bam(dir, out_name):
    """Collect information from endpoint BAM

    Inputs -
        dir - str
    Returns -
        tuple (int, int)
    """
    n_assigned_umi = 0

    with gzip.open(out_name, 'wt') as fh:
        fh.write('read_name\tcell_barcode\tumi\n')
        try:
            bam = pysam.AlignmentFile(f'{dir}/final_files/03.umis/passedParsed.bam', 'rb')
            for read in bam.fetch():
                # Count primary alignments
                if not read.is_secondary and not read.is_supplementary:
                    # Read name
                    name = get_original_name(read.query_name)

                    # Auxiliary tags
                    tags = read.get_tags()

                    bc = None
                    umi = None
                    for tag in tags:
                        id = tag[0]
                        val = tag[1]

                        if id == 'BC':
                            bc = val
                        elif id == 'U8':
                            umi = val
                            n_assigned_umi += 1

                        if bc is not None and umi is not None:
                            break

                    if bc is None:
                        bc = 'NA'
                    if umi is None:
                        umi = 'NA'

                    fh.write(f'{name}\t{bc}\t{umi}\n')

            bam.close()
        except FileNotFoundError:
            # It's okay if the file isn't found - it just means sicelore couldn't
            # generate this file, so leave n_assigned_umi as 0
            pass

    return n_assigned_umi

def main():
    # Prep output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Directories to extract information from
    DIRS = [
        '100_500bp_005M',
        '100_500bp_025M',
        '100_500bp_050M',
        '100_500bp_075M',
        '100_500bp_100M',
        '3p_repeat',
        '5p_repeat',
        'concat_2x',
        'concat_2x_rev',
        'concat_2x_rev_comp',
        'concat_3x',
        'fwd_fwd_rev',
        'fwd_rev_fwd',
        'fwd_rev_rev',
        'rev_fwd_fwd',
        'rev_fwd_rev',
        'rev_rev_fwd',
        'A375_10mil',
        '500_1000bp_10M',
        '1000_1500bp_10M',
        '1500_2000bp_10M',
        '2000_2500bp_10M',
    ]

    out = {
        'dataset': [],
        't_step_1': [],
        't_step_2': [],
        't_step_3': [],
        'max_memory_kb': [],
        'n_arch_passed': [],
        'n_arch_failed': [],
        'n_assigned_bc': [],
        'n_primary': [],
        'n_assigned_umi': [],
    }
    for dir in DIRS:
        print(f'Processing {dir}')
        df = get_pipeline_info(dir)

        t_begin = time.time()

        # Number of reads passing initial filters
        print('\tcounting passed reads ... ')
        t_start = time.time()
        n_passed, n_has_bc, n_failed = get_passed_reads(dir, df)
        t_end = time.time()
        print(f'\t... done ({t_end-t_start:.6f} sec)')

        # Timing
        print('\tcollecting time information', end=' ... ')
        t_start = time.time()
        t1, t2, t3 = get_timing(df)
        t_end = time.time()
        print(f'done ({t_end-t_start:.6f} sec)')

        # Maximum memory usage
        print('\tfinding maximum memory usage', end=' ... ')
        t_start = time.time()
        max_mem = get_peak_memory(df)
        t_end = time.time()
        print(f'done ({t_end-t_start:.6f} sec)')

        # Get number of primary alignments from minimap2 BAM
        print('\tcounting number of primary alignments', end=' ... ')
        t_start = time.time()
        n_primary = get_n_primary_alignments(dir)
        t_end = time.time()
        print(f'done ({t_end-t_start:.6f} sec)')

        # Parse BAM file
        print('\treading BAM file', end=' ... ')
        t_start = time.time()
        n_assigned_umi = process_bam(dir, f'{OUTPUT_DIR}/barcodes.{dir}.tsv.gz')
        t_end = time.time()
        print(f'done ({t_end-t_start:.6f} sec)')

        t_final = time.time()
        print(f'\tcompleted entire process in {t_final-t_begin:.6f} sec')

        out['dataset'].append(dir)
        out['t_step_1'].append(t1)
        out['t_step_2'].append(t2)
        out['t_step_3'].append(t3)
        out['max_memory_kb'].append(max_mem)
        out['n_arch_passed'].append(n_passed)
        out['n_arch_failed'].append(n_failed)
        out['n_assigned_bc'].append(n_has_bc)
        out['n_primary'].append(n_primary)
        out['n_assigned_umi'].append(n_assigned_umi)

    data = pd.DataFrame(out)

    data.to_csv(
        f'{OUTPUT_DIR}/benchmarking_results.tsv',
        sep='\t',
        na_rep='NA',
        index=False
    )

    # Write README describing benchmarking results files
    write_readme()

    return None

if __name__ == '__main__':
    main()
