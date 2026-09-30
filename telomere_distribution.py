#!/usr/bin/env python

import argparse
import re
import sys
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="Input FASTA file", required=True)
parser.add_argument("-l", "--length", type=int, help="Scan the first/last N bp (default = 300 bp)", required=False)
parser.add_argument("-f", "--forward", type=str, help="Regular expression to search start of contig (default = TTTAGGG)")
parser.add_argument("-r", "--reverse", type=str, help="Regular expression to search end of contig (default = CCCTAAA)")

args = parser.parse_args()

# Default window size (1000 bp)
window_size = 1000

# Regular expressions for telomere motifs
telomere_F = "TTTAGGG"  # Forward motif
telomere_R = "CCCTAAA"  # Reverse motif

if args.forward:
    telomere_F = args.forward
if args.reverse:
    telomere_R = args.reverse

# Process each sequence in the input FASTA file
for record in SeqIO.parse(args.input, "fasta"):
    sequence_id = str(record.description)
    sequence = str(record.seq).upper().strip("N")
    contig_length = len(sequence)

    # Initialize counters for each window
    forward_counts = []
    reverse_counts = []

    # Loop over the sequence in steps of `window_size`
    for start in range(0, contig_length, window_size):
        end = min(start + window_size, contig_length)
        window_seq = sequence[start:end]

        # Count forward and reverse telomeres in the window
        forward_count = len(re.findall(telomere_F, window_seq))
        reverse_count = len(re.findall(telomere_R, window_seq))

        forward_counts.append(forward_count)
        reverse_counts.append(reverse_count)

    # Generate a plot for this contig
    plt.figure(figsize=(12, 6))

    # Plot both forward and reverse telomere counts
    plt.plot(range(0, contig_length, window_size), forward_counts, label='Forward TTTAGGG', color='blue', marker='o')
    plt.plot(range(0, contig_length, window_size), reverse_counts, label='Reverse CCCTAAA', color='red', marker='x')

    # Add labels and title
    plt.title(f'Distribution des Télomères sur le Contig {sequence_id}')
    plt.xlabel('Position sur le Contig (fenêtres de 1000 bases)')
    plt.ylabel('Nombre de motifs télomériques')
    plt.legend()

    # Save the plot as an image
    plt.savefig(f'telomere_distribution_{sequence_id}.png')
    plt.close()  # Close the plot to avoid overlap

    # Print summary of telomere counts per contig
    print(f"Contig {sequence_id}:")
    print(f"  Forward TTTAGGG counts: {sum(forward_counts)}")
    print(f"  Reverse CCCTAAA counts: {sum(reverse_counts)}")

print("Graphs saved for all sequences.")
