#!/usr/bin/env python3

import os
import csv
import subprocess
from Bio import AlignIO
import numpy as np
import pandas as pd

# ---- TXT CONCATENATION FUNCTION ----
def convert_txts_in_va_env_to_single_csv(base_directory, output_csv_path, delimiter="\t"):
    all_data = []
    header = None
    for entry in os.listdir(base_directory):
        va_env_dir = os.path.join(base_directory, entry)
        if os.path.isdir(va_env_dir) and entry.startswith("VA_ENV"):
            for root, dirs, files in os.walk(va_env_dir):
                for file_name in files:
                    if file_name.endswith(".txt"):
                        file_path = os.path.join(root, file_name)
                        with open(file_path, 'r', encoding='utf-8') as txt_file:
                            lines = txt_file.readlines()
                            if len(lines) < 2:
                                continue
                            file_header = lines.strip().split(delimiter)
                            file_data = [line.strip().split(delimiter) for line in lines[1:]]
                            if header is None:
                                header = file_header
                            elif header != file_header:
                                raise ValueError(f"Header mismatch in file: {file_path}")
                            all_data.extend(file_data)
    if header is None or not all_data:
        print("No valid data found to write.")
        return
    with open(output_csv_path, 'w', newline='', encoding='utf-8') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(header)
        writer.writerows(all_data)
    print(f"Concatenation completed! Single CSV saved as '{output_csv_path}'.")

# ---- SEQUENCE ANALYSIS FUNCTIONS ----
def align_symbol_sequences(input_folder, output_file, symbol):
    fasta_files = [f for f in os.listdir(input_folder) if f.endswith('.fasta')]
    combined_fasta = os.path.join(input_folder, f'combined_{symbol}_sequences.fasta')
    with open(combined_fasta, 'w') as outfile:
        for fasta_file in fasta_files:
            file_path = os.path.join(input_folder, fasta_file)
            with open(file_path, 'r') as infile:
                header = ''
                seq_lines = []
                for line in infile:
                    if line.startswith('>'):
                        if seq_lines and symbol in header:
                            outfile.write(header)
                            outfile.writelines(seq_lines)
                        header = line
                        seq_lines = []
                    else:
                        seq_lines.append(line)
                if seq_lines and symbol in header:
                    outfile.write(header)
                    outfile.writelines(seq_lines)
    mafft_cmd = ['mafft', '--auto', combined_fasta]
    with open(output_file, 'w') as out:
        subprocess.run(mafft_cmd, stdout=out, check=True)
    print(f'MAFFT alignment saved to {output_file}')

def make_snp_matrix(alignment_file, output_file):
    alignment = AlignIO.read(alignment_file, "fasta")
    n = len(alignment)
    dist_matrix = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            seq1 = alignment[i].seq
            seq2 = alignment[j].seq
            snps = sum((b1 != b2) and (b1 != '-') and (b2 != '-') for b1, b2 in zip(seq1, seq2))
            dist_matrix[i, j] = dist_matrix[j, i] = snps
    labels = [record.id for record in alignment]
    df = pd.DataFrame(dist_matrix, index=labels, columns=labels)
    df.to_csv(output_file)
    print(f"SNP distance matrix written to {output_file}")

def make_phylo_tree(alignment_file, tree_file):
    if not os.path.exists(alignment_file):
        raise FileNotFoundError(f"Alignment file not found: {alignment_file}")
    try:
        with open(tree_file, 'w') as outtree:
            subprocess.run(['fasttree', '-nt', alignment_file], stdout=outtree, check=True)
        print(f"Tree saved to {tree_file}")
    except Exception as e:
        raise RuntimeError(f"FastTree failed: {str(e)}")

# ---- MAIN WORKFLOW ----
if __name__ == "__main__":
    # TXT concatenation step
    base_dir = "/home/tasmine/scratch2/AMR/results"
    output_csv = os.path.join(base_dir, "concatenated_output.csv")
    convert_txts_in_va_env_to_single_csv(base_dir, output_csv)

    # Sequence analysis step
    input_folder = '/home/tasmine/scratch2/AMR/amr_results/fast_files'
    target_symbol = "YOUR_SYMBOL"   # Set your gene or element here, e.g., 'tet(S)', 'erm(B)', etc.
    # Filenames using symbol
    aligned_fasta = f'/home/tasmine/scratch2/AMR/aligned_{target_symbol}.fasta'
    snp_matrix_file = f'/home/tasmine/scratch2/AMR/snp_matrix_{target_symbol}.csv'
    tree_file = f'/home/tasmine/scratch2/AMR/tree_{target_symbol}.nwk'
    # Run analysis
    align_symbol_sequences(input_folder, aligned_fasta, target_symbol)
    make_snp_matrix(aligned_fasta, snp_matrix_file)
    make_phylo_tree(aligned_fasta, tree_file)
