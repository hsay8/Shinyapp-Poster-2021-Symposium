import os
import sys

#Get classification from gtdbtk results summary
def get_classification(fasta_name):
    try:
        with open("./gtdbtk/gtdbtk.bac120.summary.tsv", "r") as gtdb_file:
            for line in gtdb_file:
                columns = line.strip().split("\t") # Turn each line into list, split by column
                if columns[0] == fasta_name:
                    classification = columns[1].split(";")  # Split the classifications so that each level of taxonomy is its own column.
                    for level in reversed(classification): # Search for the deepest classification, and return it
                        if len(level.strip()) > 3: #If there's only 3 characters, the classification is blank - so only return the first string that is longer than 3 characters. 
                            class1 = level.split("__")[1]
                            class2_mapping = {"s": "Species", "g": "Genus", "f": "Family", "o": "Order", "c" : "Class", "p" : "Phylum", "k" : "Kingdom"}
                            class2 = class2_mapping.get(level.split("__")[0], "No Classification")
                            return class1, class2
    except FileNotFoundError:
        print(f"File not found for {fasta_name}. Check GTDBTK output")

#Calculate the size of the assembly by its fasta by summing the number of bases 
def get_size(fasta_name):
    try:
        with open(f"{fasta_name}.fasta", "r") as fasta_file:
            content = fasta_file.read().splitlines()[1:] # Read in file, make sure to skip the first line (its metadata)
            size = sum(len(line) for line in content)
            return size
    except FileNotFoundError:
        print(f"Fasta file not found for {fasta_name}.")

#Get coverage from mosdepth results
def get_coverage(fasta_name):
    try:
        with open(f"./mosdepth/polished_depth/{fasta_name}.mosdepth.summary.txt") as mosdepth_file:
            for line in mosdepth_file.readlines()[1:-3]:
                columns = line.strip().split()
                return columns[-3]
    except FileNotFoundError:
        print(f"Coverage file not found for {fasta_name}. Check mosdepth output.")

#Get completeness from checkm results
def get_completeness(fasta_name):
    try:
        with open("./checkm/completeness.tab", "r") as checkm_file:
            for line in checkm_file.readlines()[3:-2]:
                columns = line.strip().split()
                if columns[0] == fasta_name:
                    return columns[-3]
    except FileNotFoundError:
        print(f"File not found for {fasta_name}. Check CheckM output.")
    
#Get contamination from checkm results
def get_contamination(fasta_name):
    try:
        with open("./checkm/completeness.tab", "r") as checkm_file:
            for line in checkm_file.readlines()[3:-2]:
                columns = line.strip().split()
                if columns[0] == fasta_name:
                    return columns[-2]
    except FileNotFoundError:
        print(f"File not found for {fasta_name}. Check CheckM output")

#Caulcuate GC 
def calculate_gc_percentage(fasta_name):
    try:
        with open(f"{fasta_name}.fasta", "r") as fasta_file:
            content = fasta_file.read().splitlines()[1:]
            total_bases = sum(len(line) for line in content if not line.startswith('>'))
            gc_bases = sum(line.count('G') + line.count('C') for line in content)
            return (gc_bases / total_bases) * 100 if total_bases > 0 else 0
    except FileNotFoundError:
        print(f"Fasta file not found for {fasta_name}.")

#Get number of rRNA genes from bakta results
def get_rRNA_genes(fasta_name):
    try: 
        with open(f"./bakta/{fasta_name}.txt", "r") as bakta_file:
            for line in bakta_file:
                if "rRNAs" in line:
                    return line.strip().split(":")[-1]
    except FileNotFoundError:
        print(f"Text file not found for {fasta_name}. Check Bakta output")

#Get number of tRNA genes from bakta results
def get_tRNA_genes(fasta_name):
    try:
        with open(f"./bakta/{fasta_name}.txt", "r") as bakta_file:
            for line in bakta_file:
                if "tRNAs" in line:
                    return line.strip().split(":")[-1]
    except FileNotFoundError:
        print(f"Text file not found for {fasta_name}. Check Bakta output")

def get_sample_name(fasta_name):
    sample = fasta_name.split("+")[0]
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
    return sample_mapping.get(sample.strip(), "No Sample")

def get_edge_num(fasta_name):
    return fasta_name.split("+")[1]

with open ("./circos_plotting/mag_summary_stats.tab", "w") as output_file:
    # Get FASTA files in the directory (sans file extension)
    fasta_files = [f[:-6] for f in os.listdir() if f.endswith(".fasta")]

    sys.stdout = output_file

    # Print header
    print("Sample\tEdge_num\tClassification\tSize\tCoverage\tCompleteness\tContamination\tGC\ttRNA_genes\trRNA_genes\tClass2\tFilename")

    for fasta_name in fasta_files:
        # Retrieve information for each column
        sample = get_sample_name(fasta_name)
        edge = get_edge_num(fasta_name)
        classification, class2 = get_classification(fasta_name)
        size = get_size(fasta_name)
        coverage = get_coverage(fasta_name)
        completeness = get_completeness(fasta_name)
        contamination = get_contamination(fasta_name)
        gc_percentage = calculate_gc_percentage(fasta_name)
        tRNA_genes = get_tRNA_genes(fasta_name)
        rRNA_genes = get_rRNA_genes(fasta_name)
        filename = fasta_name

        print(f"{sample}\t{edge}\t{classification}\t{size}\t{coverage}\t{completeness}\t{contamination}\t{gc_percentage:.2f}\t{tRNA_genes}\t{rRNA_genes}\t{class2}\t{filename}")
