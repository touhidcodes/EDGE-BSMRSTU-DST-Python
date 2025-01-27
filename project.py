import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from fpdf import FPDF

# Function to load a FASTA file and return sequences as a dictionary
def getSequences(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def main():
    file_path = input("Enter the path to your FASTA file: ").strip()

    if not os.path.isfile(file_path):
        print("File not found. Please check the file path and try again.")
        return
    
    # Load sequences from the FASTA file
    sequences = getSequences(file_path)
    if not sequences:
        print("No sequences found in the FASTA file.")
        return

    print(sequences)

main()