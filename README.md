SARS-CoV-2 Genome Analysis

Introduction

This repository focuses on bioinformatics analysis of the SARS-CoV-2 genome. The goal is to identify protein-coding genes, validate predictions using well-characterized genomes, and analyze functional characteristics such as hydrophobicity to predict transmembrane proteins.

Setup Instructions

1. Set Up Environment

Create a working directory:

mkdir sars_cov2_analysis
cd sars_cov2_analysis

2. Install Required Packages

Ensure Python and necessary dependencies are available:

module load python
pip install biopython matplotlib jupyterlab seaborn pandas

Gene Prediction and Validation

Download and Analyze the SARS-CoV-2 Genome

Fetch the SARS-CoV-2 genome using Biopython:

from Bio import Entrez
Entrez.email = "your.email@charlotte.edu"
handle = Entrez.efetch(db="nucleotide", id="NC_045512", rettype="fasta", retmode="text")
with open("sars_cov2.fasta", "w") as output_file:
    output_file.write(handle.read())
handle.close()

Identify Open Reading Frames (ORFs)

Use the following function to identify ORFs in the SARS-CoV-2 genome:

from Bio import SeqIO

def find_orfs(seq_record, min_length=30):
    """Identify ORFs in a given sequence."""
    ... # ORF detection code here

Protein Hydrophobicity Analysis

Generate a Protein FASTA File

Export ORFs as protein sequences:

def export_trimmed_proteins(orfs, output_filename):
    """Save protein sequences in FASTA format."""
    ... # Protein export function here

Compute Hydrophobicity

Load hydrophobicity values and analyze protein sequences:

import pandas as pd
hydro_table = load_hydrophobicity_table("aminoacid_properties.csv")
sars_cov2_data = calculate_hydrophobicities("sars_cov2_proteome.fasta", hydro_table)
df = pd.DataFrame(sars_cov2_data)

Comparative Viral Genome Analysis

Download Viral Metadata

Visit NCBI Virus

Download a randomized subset (20 records per host category) as a CSV file

Rename it to viral_metadata.csv and upload it to the working directory

Generate a Comparative Plot

python viral_genome_histogram.py

License

This project is licensed under the CC BY-NC-SA 4.0 License.
