from pathlib import Path
import gzip
from collections import defaultdict
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import os

# ====== CONFIG ======
PROJECT = Path(__file__).resolve().parent
FASTQ_PATH = os.path.join(PROJECT, "data", "SRR36584168.fastq")  
MOTIF = "GCGATCAACTCGCTGG"
N_FIRST = 1000
MAX_READS_FOR_PLOTS = 200000  
# =====================

def open_maybe_gzip(path):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8", errors="replace")

def count_overlapping(subseq: str, seq: str) -> int:
    count = 0
    start = 0
    while True:
        idx = seq.find(subseq, start)
        if idx == -1:
            return count
        count += 1
        start = idx + 1

def main():
    read_count = 0
    first_id = first_seq = first_qual = None
    motif_total = 0

    out_first1000 = os.path.join(PROJECT, "data", "first_1000.fastq.gz")
    out_handle = gzip.open(out_first1000, "wt")

    per_pos_quals = defaultdict(list)
    lengths = []

    with open_maybe_gzip(FASTQ_PATH) as handle:
        for rec in SeqIO.parse(handle, "fastq"):
            read_count += 1
            seq = str(rec.seq)
            quals = rec.letter_annotations["phred_quality"]

            if read_count == 1:
                first_id = rec.id
                first_seq = seq
                first_qual = quals

            motif_total += count_overlapping(MOTIF, seq)

            if read_count <= N_FIRST:
                SeqIO.write(rec, out_handle, "fastq")

            if read_count <= MAX_READS_FOR_PLOTS:
                lengths.append(len(seq))
                for i, q in enumerate(quals):
                    per_pos_quals[i].append(q)

    out_handle.close()

    # I–IV
    print("I) Number of reads:", read_count)
    print("\nII) First read")
    print("Identifier:", first_id)
    print("Sequence:", first_seq)
    print("Quality (phred):", first_qual)
    print(f"\nIII) Total occurrences of motif {MOTIF} (overlapping counted):", motif_total)
    print("\nIV) Wrote first 1000 reads to:", out_first1000)

    # V: per-base quality boxplot
    if per_pos_quals:
        
        max_pos = max(per_pos_quals.keys()) if per_pos_quals else -1
        data = [per_pos_quals[i] for i in range(max_pos + 1)]

        plt.figure(figsize=(16, 6))
        plt.boxplot(data, showfliers=False)

        plt.xlabel("Base position (1-based)")
        plt.ylabel("Phred quality")
        plt.title("Per-base quality (boxplot)")

        # Show fewer x tick labels to avoid overlap
        step = 10 if (max_pos + 1) <= 200 else 25
        tick_positions = np.arange(1, max_pos + 2, step)       
        tick_labels = [str(x) for x in tick_positions]
        plt.xticks(tick_positions, tick_labels, rotation=0)

        plt.grid(True, axis="y", alpha=0.2)

        quality_plot = os.path.join(PROJECT, "plots", "per_base_quality_boxplot.png")
        plt.savefig(quality_plot, dpi=200, bbox_inches="tight")
        plt.close()
        print("\nV) Saved quality boxplot:", quality_plot)


    # VI: read length density plot
    if lengths:
       
        lengths_arr = np.array(lengths)

        plt.figure(figsize=(12, 6))
        plt.hist(lengths_arr, bins=60, density=True)

        plt.xlabel("Read length")
        plt.ylabel("Density")
        plt.title("Read length distribution (density)")

        # Visual enhancement: 
        plt.yscale("log")

        plt.grid(True, alpha=0.2)

        len_plot = os.path.join(PROJECT, "plots", "read_length_density.png")
        plt.savefig(len_plot, dpi=200, bbox_inches="tight")
        plt.close()
        print("VI) Saved read length density plot:", len_plot)


if __name__ == "__main__":
    main()
