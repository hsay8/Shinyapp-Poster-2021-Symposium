import os
import sys

def get_edge_name(fasta_name):
    try:
        with open(f"{fasta_name}.fasta", "r") as fasta_file:
            edge = fasta_file.readline().split(">")[-1] # Read in file first line, but not including the > symbol
            return edge.strip()
    except FileNotFoundError:
        print(f"Fasta file not found for {fasta_name}.")

def get_contig_names(fasta_name):
    sample_mapping = {  "8-carbon_granules4": "GAC0",
                        "7-T2_jan27-20" : "GAC1",
                        "3-T4_feb14-20": "GAC2",
                        "9-gac_march_02-20": "GAC3",
                        "6-T2_oct6-20": "GAC4",
                        "5-T2_may3-21": "GAC5",
                        "4-T3_may7-21": "GAC6",
                        "2-T4_may14-21": "GAC7",
                        "1-T1_aug10-21": "GAC8",
                        "10-T3_apr21-22": "GAC9"     }
    
    for key, value in sample_mapping.items():
        if key in fasta_name:
            return fasta_name.replace(key, value)
    return fasta_name

def get_size(fasta_name):
    try:
        with open(f"{fasta_name}.fasta", "r") as fasta_file:
            content = fasta_file.read().splitlines()[1:] # Read in file, make sure to skip the first line (its metadata)
            size = sum(len(line) for line in content)
            return size
    except FileNotFoundError:
        print(f"Fasta file not found for {fasta_name}.")

with open("cytoband.txt", "w") as output_file:
    fasta_files = [f[:-6] for f in os.listdir() if f.endswith(".fasta")]

    sys.stdout = output_file

    for fasta_name in fasta_files:
        edge = get_edge_name(fasta_name)
        name = get_contig_names(fasta_name)
        size = get_size(fasta_name)
        print(f"{edge}\t{size}\tnot real\tnot real\t{name}")

