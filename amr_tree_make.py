#!/usr/bin/env python3

import os
import subprocess
from Bio import AlignIO
import numpy as np
import pandas as pd

# ---- USER SETTINGS: Change paths as needed ----
input_folder = '/home/tasmine/scratch2/AMR/amr_results'        # folder with input FASTA files
target_symbol = "YOUR_SYMBOL"   # e.g., 'tet(S)', 'erm(B)', etc.
aligned_fasta = f'/home/tasmine/scratch2/AMR/aligned_{target_symbol}.fasta'
snp_matrix_file = f'/home/tasmine/scratch2/AMR/snp_matrix_{target_symbol}.csv'
tree_file = f'/home/tasmine/scratch2/AMR/tree_{target_symbol}.nwk'


# Dynamically name files using the target symbol
aligned_fasta = f'/home/tasmine/scratch2/AMR/aligned_{target_symbol}.fasta'
snp_matrix_file = f'/home/tasmine/scratch2/AMR/snp_matrix_{target_symbol}.csv'
tree_file = f'/home/tasmine/scratch2/AMR/tree_{target_symbol}.nwk'

def align_symbol_sequences(input_folder, output_file, symbol):
    """Extract sequences for a given symbol, combine, align with MAFFT."""
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
                        # Write previous record if it's the target symbol
                        if seq_lines and symbol in header:
                            outfile.write(header)
                            outfile.writelines(seq_lines)
                        header = line
                        seq_lines = []
                    else:
                        seq_lines.append(line)
                # Write last record of the file
                if seq_lines and symbol in header:
                    outfile.write(header)
                    outfile.writelines(seq_lines)

    # Run MAFFT to align the extracted sequences
    mafft_cmd = ['mafft', '--auto', combined_fasta]
    with open(output_file, 'w') as out:
        subprocess.run(mafft_cmd, stdout=out, check=True)
    print(f'MAFFT alignment saved to {output_file}')

def make_snp_matrix(alignment_file, output_file):
    """Create SNP distance matrix from aligned FASTA."""
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
    """Construct a phylogenetic tree using FastTree."""
    if not os.path.exists(alignment_file):
        raise FileNotFoundError(f"Alignment file not found: {alignment_file}")
    try:
        with open(tree_file, 'w') as outtree:
            subprocess.run(['fasttree', '-nt', alignment_file], stdout=outtree, check=True)
        print(f"Tree saved to {tree_file}")
    except Exception as e:
        raise RuntimeError(f"FastTree failed: {str(e)}")

# ------ MAIN WORKFLOW -------
if __name__ == "__main__":
    align_symbol_sequences(input_folder, aligned_fasta, target_symbol)
    make_snp_matrix(aligned_fasta, snp_matrix_file)
    make_phylo_tree(aligned_fasta, tree_file)
