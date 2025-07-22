import os
import re
import csv
import polars as pl
import math
import typer
import logging
import random
import numpy as np
from typing import List
from multiprocessing import Pool
from Bio import SeqIO
from Levenshtein import distance as lev_distance

########################## Sequence Order Loader ##########################

def seq_orders(file_path, model):
    try:
        if not os.path.isfile(file_path):
            print(f"The file '{file_path}' does not exist.")
            return [], [], [], []

        with open(file_path, 'r') as file:
            for line in file:
                fields = line.strip().split("\t")
                if len(fields) < 5:
                    continue
                model_name = fields[0].strip()
                sequence_order = [seg.strip() for seg in fields[1][1:-1].split(',')]
                sequences = [seg.strip() for seg in fields[2][1:-1].split(',')]
                barcodes = fields[3].strip().strip("'").strip('"').split(',')
                UMIs = fields[4].strip().strip("'").strip('"').split(',')
                if model_name == model:
                    return sequence_order, sequences, barcodes, UMIs
    except Exception as e:
        print(f"An error occurred: {e}")
    return [], [], [], []

########################## Segment Generation ##########################

def mutate_sequence(seq: str, num_mutations: int) -> str:
    seq = list(seq)
    positions = random.sample(range(len(seq)), num_mutations)
    for pos in positions:
        seq[pos] = random.choice([b for b in "ACGT" if b != seq[pos]])
    return ''.join(seq)

def mutate_to_target_ld(original: str, whitelist: list, target_ld: int, max_attempts: int = 30) -> str:
    """Mutate a given CBC to have a specific LD to its nearest whitelist entry."""
    if target_ld == 0:
        return original

    for _ in range(max_attempts):
        mutated = mutate_sequence(original, target_ld)
        if mutated in whitelist:
            continue
        min_dist = min(lev_distance(mutated, w) for w in whitelist)
        if min_dist == target_ld:
            return mutated
    return original  # fallback to original if unable

def introduce_errors_with_labels_context(sequence, label, mismatch_rate, insertion_rate, deletion_rate,
                                         cbc_error_rate, polyT_error_rate, max_insertions, whitelist=None):
    error_sequence, error_labels = [], []
    i = 0

    # Desired LD histogram for CBCs (peak at 0)
    ld_distribution = {
        0: 0.45,
        1: 0.25,
        2: 0.15,
        3: 0.08,
        4: 0.05,
        5: 0.02
    }
    ld_vals, ld_probs = zip(*ld_distribution.items())

    while i < len(sequence):
        base = sequence[i]
        lbl = label[i]
        insertion_count = 0

        # Handle CBC region
        if lbl == "CBC" and whitelist is not None:
            cbc_seq = []
            start = i
            while i < len(sequence) and label[i] == "CBC":
                cbc_seq.append(sequence[i])
                i += 1
            cbc_seq = ''.join(cbc_seq)
            target_ld = np.random.choice(ld_vals, p=ld_probs)
            mutated_cbc = mutate_to_target_ld(cbc_seq, whitelist, target_ld)
            error_sequence.extend(mutated_cbc)
            error_labels.extend(["CBC"] * len(mutated_cbc))
            continue

        # General logic for other regions
        r = np.random.random()
        rates = (polyT_error_rate if lbl in ("polyT", "polyA") else mismatch_rate,
                 polyT_error_rate if lbl in ("polyT", "polyA") else insertion_rate,
                 polyT_error_rate if lbl in ("polyT", "polyA") else deletion_rate)

        if lbl == 'ACC':
            rates = (0, 0, 0)
        elif lbl == 'CBC':
            rates = (cbc_error_rate, cbc_error_rate, cbc_error_rate)

        if r < rates[0]:
            error_sequence.append(np.random.choice([b for b in "ATCG" if b != base]))
            error_labels.append(lbl)
        elif r < rates[0] + rates[1]:
            error_sequence.append(base)
            error_labels.append(lbl)
            while insertion_count < max_insertions:
                error_sequence.append(np.random.choice(list("ATCG")))
                error_labels.append(lbl)
                insertion_count += 1
                if np.random.random() >= rates[1]:
                    break
        elif r < sum(rates):
            pass  # deletion
        else:
            error_sequence.append(base)
            error_labels.append(lbl)

        i += 1

    return "".join(error_sequence), error_labels

def generate_segment(segment_type, segment_pattern, length_range, transcriptome_records=None, cbc_pool=None):
    if segment_type.lower() == "cbc" and cbc_pool:
        sequence = random.choice(cbc_pool)
        label = [segment_type] * len(sequence)
    elif re.match(r"N\d+", segment_pattern):
        length = int(segment_pattern[1:])
        sequence = "".join(np.random.choice(list("ATCG")) for _ in range(length))
        label = [segment_type] * length
    # elif segment_pattern in ["NN", "RN"] and segment_type == "cDNA" and transcriptome_records:
    #     length = np.random.randint(length_range[0], length_range[1])
    #     transcript = random.choice(transcriptome_records)
    #     transcript_seq = str(transcript.seq)
    #     if len(transcript_seq) > length:
    #         fragment = transcript_seq[:length] if random.random() < 0.5 else transcript_seq[-length:]
    #     else:
    #         fragment = transcript_seq
    #     sequence = fragment
    #     label = ["cDNA"] * len(sequence)
    elif segment_pattern in ["NN", "RN"] and segment_type == "cDNA" and transcriptome_records:
        length = np.random.randint(length_range[0], length_range[1])
        while True:
            transcript = random.choice(transcriptome_records)
            transcript_seq = str(transcript.seq)
            if len(transcript_seq) >= length:
                break
        fragment = transcript_seq[:length] if random.random() < 0.5 else transcript_seq[-length:]
        sequence = fragment
        label = ["cDNA"] * len(sequence)
    elif segment_pattern in ["A", "T"]:
        length = np.random.randint(0, 50)
        sequence = segment_pattern * length
        label = ["polyA" if segment_pattern == "A" else "polyT"] * length
    else:
        sequence = segment_pattern
        label = [segment_type] * len(sequence)
    return sequence, label

########################## Utility Functions ##########################

def reverse_complement(sequence):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement[base] for base in reversed(sequence))

def reverse_labels(labels):
    return labels[::-1]

def simulate_quality_scores(read_len, platform="ONT"):
    platform_q_params = {
        "ONT": (20, 4),
        "PacBio": (20, 3),
        "Illumina": (30, 5),
        "BGI": (28, 4),
        "IonTorrent": (25, 6)
    }
    mean_q, std_q = platform_q_params.get(platform, (30, 5))
    return np.clip(np.random.normal(loc=mean_q, scale=std_q, size=read_len).round(), 2, 40).astype(int).tolist()

def generate_unique_barcodes(n, k=16, min_dist=3):
    barcodes = []
    tries = 0
    while len(barcodes) < n and tries < n * 100:
        candidate = "".join(random.choices("ACGT", k=k))
        if all(lev_distance(candidate, b) >= min_dist for b in barcodes):
            barcodes.append(candidate)
        tries += 1
    return barcodes

def reverse_or_revcomp_sequence_list(seq_list, operation):
    """Reverse or reverse-complement a list of sequences."""
    transformed = []
    for s in seq_list[::-1]:
        if re.match(r"N\d+", s):  
            transformed.append(s)
        elif operation == "revcomp":
            transformed.append(reverse_complement(s))
        elif operation == "reverse":
            transformed.append(s[::-1])
        else:
            transformed.append(s)
    return transformed

########################### Chunked Processing ##########################

def simulate_chunk_reads(chunk_id, num_reads, mismatch_rate, insertion_rate, deletion_rate, cbc_error_rate,
                         polyT_error_rate, max_insertions, segments_order, segments_patterns,
                         length_range, rc, transcriptome_records, cbc_pool, output_dir, platform, max_extra_bases):
    reads, labels, metadata = [], [], []
    max_reads_per_file = 4000

    for read_id in range(num_reads):
        full_id = f"chunk{chunk_id}_read{read_id}"
        read_segments, label_segments, read_meta = [], [], {"read_id": full_id}

        for seg_type, seg_pattern in zip(segments_order, segments_patterns):
            segment_seq, segment_label = generate_segment(seg_type, seg_pattern, length_range, transcriptome_records, cbc_pool)
            read_segments.append(segment_seq)
            label_segments.append(segment_label)
            if seg_type.lower() in ["cbc", "i5", "i7", "umi"]:
                read_meta[seg_type.upper()] = segment_seq

        # Add extra bases to start/end
        extra_start = "".join(random.choices("ACGT", k=random.randint(0, max_extra_bases)))
        extra_end = "".join(random.choices("ACGT", k=random.randint(0, max_extra_bases)))
        error_read_segments, error_label_segments = [], []

        for segment, segment_label in zip(read_segments, label_segments):
            error_segment, error_labels = introduce_errors_with_labels_context(
                segment, segment_label, mismatch_rate, insertion_rate, deletion_rate, cbc_error_rate, polyT_error_rate, max_insertions)
            error_read_segments.append(error_segment)
            error_label_segments.extend(error_labels)

        read = extra_start + "".join(error_read_segments) + extra_end
        reads.append(read)
        labels.append(["extra"] * len(extra_start) + error_label_segments + ["extra"] * len(extra_end))
        metadata.append(read_meta)

        if rc:
            rc_read = reverse_complement(read)
            rc_labels = reverse_labels(labels[-1])
            reads.append(rc_read)
            labels.append(rc_labels)
            metadata.append({**read_meta, "read_id": f"{full_id}_rc"})

    os.makedirs(os.path.join(output_dir, "simulated_data"), exist_ok=True)
    total_reads = len(reads)
    num_files = math.ceil(total_reads / max_reads_per_file)

    for i in range(num_files):
        fq_path = os.path.join(output_dir, f"simulated_data/reads_chunk{chunk_id}_{i}.fastq")
        with open(fq_path, 'w') as fq:
            start = i * max_reads_per_file
            end = min((i + 1) * max_reads_per_file, total_reads)
            for j in range(start, end):
                q = simulate_quality_scores(len(reads[j]), platform)
                fq.write(f"@{metadata[j]['read_id']}\n{reads[j]}\n+\n{''.join(chr(33 + x) for x in q)}\n")

    meta_path = os.path.join(output_dir, f"simulated_data/metadata_chunk{chunk_id}.tsv")
    with open(meta_path, 'w', newline='') as mf:
        writer = csv.DictWriter(mf, fieldnames=sorted(metadata[0].keys()), delimiter='\t')
        writer.writeheader()
        for meta in metadata:
            writer.writerow(meta)



########################### Typer App ##########################

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
app = typer.Typer()

@app.command()
def simulate_data(model_name: str,
                  output_dir: str,
                  transcriptome_fasta: str, 
                  invalid_type: str = typer.Option(None, help="what type of invalid reads to generate?"),
                  num_reads: int = typer.Option(50000, help="Number of reads to be simulated (Same reads will be reverse complemented if rc == True)"),
                  mismatch_rate: float = typer.Option(0.05, help="mismatch rate"),
                  insertion_rate: float = typer.Option(0.05, help="insertion rate"), 
                  deletion_rate: float = typer.Option(0.0612981959469103, help="deletion rate"),
                  min_cDNA: int = typer.Option(100, help="minimum cDNA length"),
                  max_cDNA: int = typer.Option(499, help="maximum cDNA length"),
                  cbc_error_rate: int = typer.Option(0, help="CBC error rate"),
                  polyT_error_rate: float = typer.Option(0.02, help="error rate within polyT or polyA segments"),
                  max_insertions: float = typer.Option(2, help="maximum number of allowed insertions after a base"),
                  threads: int = typer.Option(2, help="number of CPU threads"), 
                  platform: str = typer.Option("ONT", help="Platform"),
                  cbc_whitelist_domain: str = typer.Option(None, help="cell barcodes whitelist file to generate the pool from"),
                  num_cbcs: int = typer.Option(6500, help="Number of cell barcodes to generate a pool to select from"),
                  rc: bool = typer.Option(True, help="whether to include reverse complements of the reads in the training data.\nFinal dataset will contain twice the number of user-specified reads"),
                  reads_per_batch: int = typer.Option(5000000, help="number of reads in each batch"),
                  cbc_min_dist: int = typer.Option(2, help="If cell barcodes are randomly generated, the minimum distance each have from every other"),
                  max_extra_bases: int = typer.Option(70, help="maximum number of bases to be added on either end of the reads")):

    logger.info("Loading transcriptome fasta file")
    transcriptome_records = list(SeqIO.parse(transcriptome_fasta, "fasta"))
    logger.info("Transcriptome fasta loaded")

    os.makedirs(f"{output_dir}/simulated_data", exist_ok=True)
    
    if cbc_whitelist_domain:
        logger.info("using the whitelist domain file to generate the whitelist pool")
        whitelist_domain = pl.read_csv(cbc_whitelist_domain, has_header = False)
        whitelist_domain_list = whitelist_domain['column_1'].to_list()
        cbc_pool = random.sample(whitelist_domain_list, num_cbcs)
    else:
        logger.info("generating whitelist pool by random base generation")
        cbc_pool = generate_unique_barcodes(num_cbcs, k=16, min_dist=cbc_min_dist)
    
    with open(f"{output_dir}/simulated_data/cbc_whitelist.txt", 'w') as f:
        for cbc in cbc_pool:
            f.write(f"{cbc}\n")

    length_range = (min_cDNA, max_cDNA)
    seq_order, sequences, barcodes, UMIs = seq_orders("training_seq_orders.tsv", model_name)

    if invalid_type == "concat_2x":
        seq_order = seq_order * 2
        sequences = sequences * 2
    elif invalid_type == "concat_3x":
        seq_order = seq_order * 3
        sequences = sequences * 3
    elif invalid_type == "concat_2x_rev":
        seq_order = seq_order + seq_order[::-1]
        sequences = sequences + reverse_or_revcomp_sequence_list(sequences, "reverse")
    elif invalid_type == "concat_2x_revcomp":
        seq_order = seq_order + seq_order[::-1]
        sequences = sequences + reverse_or_revcomp_sequence_list(sequences, "revcomp")

    elif invalid_type == "fwd_rev_fwd":
        seq_order = seq_order + seq_order[::-1] + seq_order
        sequences = sequences + reverse_or_revcomp_sequence_list(sequences, "revcomp") + sequences
    elif invalid_type == "fwd_rev_rev":
        seq_order = seq_order + seq_order[::-1] + seq_order[::-1] 
        sequences = sequences + reverse_or_revcomp_sequence_list(sequences, "revcomp") + reverse_or_revcomp_sequence_list(sequences, "revcomp")
    elif invalid_type == "fwd_fwd_rev":
        seq_order = seq_order + seq_order + seq_order[::-1]
        sequences = sequences + sequences + reverse_or_revcomp_sequence_list(sequences, "revcomp")
    elif invalid_type == "rev_rev_fwd":
        seq_order = seq_order[::-1] + seq_order[::-1] + seq_order
        sequences = reverse_or_revcomp_sequence_list(sequences, "revcomp") + reverse_or_revcomp_sequence_list(sequences, "revcomp") + sequences
    elif invalid_type == "rev_fwd_fwd":
        seq_order = seq_order[::-1] + seq_order + seq_order
        sequences = reverse_or_revcomp_sequence_list(sequences, "revcomp") + sequences + sequences
    elif invalid_type == "rev_fwd_rev":
        seq_order = seq_order[::-1] + seq_order + seq_order[::-1]
        sequences = reverse_or_revcomp_sequence_list(sequences, "revcomp") + sequences + reverse_or_revcomp_sequence_list(sequences, "revcomp")
    elif invalid_type == "5p_repeat":
        seq_order = seq_order[0] * random.choice(range(3,15))
        sequences = sequences[0] * random.choice(range(3,15))
    elif invalid_type == "3p_repeat":
        seq_order = seq_order[-1] * random.choice(range(3,15))
        sequences = sequences[-1] * random.choice(range(3,15))

    batch_id = 0
    reads_generated = 0

    while reads_generated < num_reads:
        current_batch = min(reads_per_batch, num_reads - reads_generated)
        chunk_size = current_batch // threads
        remainder = current_batch % threads

        chunk_args = []
        for i in range(threads):
            n = chunk_size + (1 if i < remainder else 0)
            if n > 0:
                chunk_args.append((batch_id * threads + i, n, mismatch_rate, insertion_rate, deletion_rate,
                                   cbc_error_rate, polyT_error_rate, max_insertions, seq_order, sequences,
                                   length_range, rc, transcriptome_records, cbc_pool, output_dir, platform, max_extra_bases))

        logger.info(f"Starting batch {batch_id + 1}, total reads so far: {reads_generated}")
        with Pool(threads) as pool:
            pool.starmap(simulate_chunk_reads, chunk_args)

        reads_generated += current_batch
        batch_id += 1

    logger.info("Combining metadata files")
    with open(f"{output_dir}/simulated_data/metadata.tsv", 'w', newline='') as out_meta:
        writer = None
        for fname in sorted(os.listdir(f"{output_dir}/simulated_data")):
            path = os.path.join(output_dir, "simulated_data", fname)
            if fname.startswith("metadata_chunk") and fname.endswith(".tsv"):
                with open(path) as f:
                    reader = csv.DictReader(f, delimiter='\t')
                    if writer is None:
                        writer = csv.DictWriter(out_meta, fieldnames=reader.fieldnames, delimiter='\t')
                        writer.writeheader()
                    for row in reader:
                        writer.writerow(row)
                os.remove(path)

    logger.info("Metadata combined. Simulation complete.")

if __name__ == "__main__":
    app()