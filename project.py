import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from fpdf import FPDF


def select_fasta_file():
    # Get all FASTA files in the current directory
    fasta_files = [file for file in os.listdir() if file.endswith(".fasta")]

    # Check if no FASTA files are found
    if not fasta_files:
        print("No FASTA files found in the current directory.")
        return
     
    # Display available FASTA files
    print("\nAvailable FASTA files:")
    for index, file_name in enumerate(fasta_files, 1):
        print(f"{index}. {file_name}")
    
      # Prompt user to select a file by number
    user_choice = input("\nEnter the number of the FASTA file you want to use: ")

    # Validate user input
    if user_choice.isdigit():
        file_index = int(user_choice) - 1
        if 0 <= file_index and file_index < len(fasta_files):
            selected_file = fasta_files[file_index]
            print(f"\nYou selected: {selected_file}")
            return selected_file
        else:
            print("\nInvalid choice. Please select a valid number from the list.")
            return
    else:
        print("\nInvalid input. Please enter a number corresponding to the file.")
        return

# Function to load a FASTA file and return sequences as a dictionary
def getSequences(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def main():
    
    selected_file = select_fasta_file()
    
    # Load sequences from the FASTA file
    sequences = getSequences(selected_file)
    if not sequences:
        print("No sequences found in the FASTA file.")
        return

    print(sequences)

main()