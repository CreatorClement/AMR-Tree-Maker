#!/usr/bin/env python3

import os
import subprocess
import pandas as pd
import numpy as np
from Bio import AlignIO

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

if __name__ == "__main__":
    # Read the concatenated CSV to get unique element symbols
    base_dir = "/home/tasmine/scratch2/AMR/results"
    output_csv = os.path.join(base_dir, "concatenated_output.csv")
    input_folder = '/home/tasmine/scratch2/AMR/amr_results/fasta_files'

    df = pd.read_csv(output_csv)

    # Get unique values in the 'Element symbol' column (assumes exact column name)
    unique_symbols = df['Element symbol'].dropna().unique()

    for symbol in unique_symbols:
        # Sanitize symbol for filenames
        safe_symbol = str(symbol).replace('(', '_').replace(')', '_').replace('/', '_')
        aligned_fasta = f'/home/tasmine/scratch2/AMR/aligned_{safe_symbol}.fasta'
        snp_matrix_file = f'/home/tasmine/scratch2/AMR/snp_matrix_{safe_symbol}.csv'
        tree_file = f'/home/tasmine/scratch2/AMR/tree_{safe_symbol}.nwk'

        print(f"\nProcessing element: {symbol}")

        align_symbol_sequences(input_folder, aligned_fasta, symbol)
        make_snp_matrix(aligned_fasta, snp_matrix_file)
        make_phylo_tree(aligned_fasta, tree_file)
