# # Import necessary libraries
# import os
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.SeqUtils import gc_fraction

# # Function to load sequences from a FASTA file
# def load_fasta(file_name):
#     """
#     Load sequences from a FASTA file and return as a dictionary.
#     """
#     sequences = {}
#     for record in SeqIO.parse(file_name, "fasta"):
#         sequences[record.id] = str(record.seq)
#     return sequences

# # Function to calculate GC content of a sequence
# def calculate_gc_content(sequence):
#     """
#     Calculate GC content as a percentage using BioPython.
#     """
#     return gc_fraction(sequence) * 100

# # Function to get the reverse complement of a DNA sequence
# def reverse_complement(sequence):
#     """
#     Generate the reverse complement of a DNA sequence.
#     """
#     seq_obj = Seq(sequence)
#     return str(seq_obj.reverse_complement())

# # Function to transcribe a DNA sequence to RNA
# def transcribe_to_rna(sequence):
#     """
#     Transcribe a DNA sequence to RNA (T -> U).
#     """
#     seq_obj = Seq(sequence)
#     return str(seq_obj.transcribe())

# # Function to analyze the sequences and return a DataFrame
# def analyze_sequences(sequences):
#     """
#     Analyze sequences and store details like length, GC content, reverse complement, and RNA transcription.
#     """
#     results = []
#     for seq_id, sequence in sequences.items():
#         gc_content = calculate_gc_content(sequence)
#         rev_complement = reverse_complement(sequence)
#         transcription = transcribe_to_rna(sequence)

#         # Append analysis results to the list
#         results.append({
#             'Sequence ID': seq_id,
#             'Length': len(sequence),
#             'GC Content (%)': round(gc_content, 2),
#             'Reverse Complement (first 50)': rev_complement[:50],
#             'RNA Transcription (first 50)': transcription[:50],
#         })

#     # Convert the list to a Pandas DataFrame
#     return pd.DataFrame(results)

# # Function to plot GC content distribution
# def plot_gc_content_distribution(df):
#     """
#     Plot the distribution of GC content using Matplotlib.
#     """
#     plt.figure(figsize=(8, 6))
#     plt.hist(df['GC Content (%)'], bins=10, color='blue', edgecolor='black', alpha=0.7)
#     plt.title('GC Content Distribution')
#     plt.xlabel('GC Content (%)')
#     plt.ylabel('Frequency')
#     plt.grid(axis='y', linestyle='--', alpha=0.7)
#     plt.show()

# # Function to plot sequence length distribution
# def plot_sequence_length_distribution(df):
#     """
#     Plot the distribution of sequence lengths using Matplotlib.
#     """
#     plt.figure(figsize=(8, 6))
#     plt.hist(df['Length'], bins=10, color='green', edgecolor='black', alpha=0.7)
#     plt.title('Sequence Length Distribution')
#     plt.xlabel('Length (bp)')
#     plt.ylabel('Frequency')
#     plt.grid(axis='y', linestyle='--', alpha=0.7)
#     plt.show()

# # Main function to tie everything together
# def main():
#     # Check for FASTA files in the current directory
#     fasta_files = [f for f in os.listdir() if f.endswith(".fasta")]
#     if not fasta_files:
#         print("No FASTA files found in the current directory.")
#         return

#     # List available FASTA files for user selection
#     print("Available FASTA files:")
#     for i, file_name in enumerate(fasta_files, 1):
#         print(f"{i}. {file_name}")

#     # Let the user choose a file
#     choice = int(input("\nSelect a FASTA file by number: ")) - 1
#     if choice < 0 or choice >= len(fasta_files):
#         print("Invalid choice. Exiting.")
#         return
#     file_name = fasta_files[choice]

#     # Load the selected FASTA file
#     print(f"\nLoading sequences from {file_name}...")
#     sequences = load_fasta(file_name)

#     if not sequences:
#         print("No sequences found in the selected FASTA file.")
#         return

#     # Analyze sequences and display results
#     print("\nAnalyzing sequences...")
#     df = analyze_sequences(sequences)

#     # Show analysis results
#     print("\n--- Sequence Analysis Results ---")
#     print(df)

#     # Plot GC content and sequence length distributions
#     print("\nGenerating plots...")
#     plot_gc_content_distribution(df)
#     plot_sequence_length_distribution(df)

#     # Save analysis results to a CSV file
#     output_file = file_name.replace(".fasta", "_analysis.csv")
#     df.to_csv(output_file, index=False)
#     print(f"\nAnalysis results saved to {output_file}")

# # if __name__ == "__main__":
# #     main()
    
# main()

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from fpdf import FPDF

# Function to load a FASTA file and return sequences as a dictionary
def load_fasta(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
        print(sequences)
    return sequences

# Function to calculate GC content using BioPython
def calculate_gc_content(sequence):
    return gc_fraction(sequence) * 100

# Function to get the reverse complement of a DNA sequence using BioPython
def reverse_complement(sequence):
    seq_obj = Seq(sequence)
    return str(seq_obj.reverse_complement())

# Function to transcribe DNA to RNA (T -> U)
def transcribe_to_rna(sequence):
    seq_obj = Seq(sequence)
    return str(seq_obj.transcribe())

# Function to calculate nucleotide composition
def nucleotide_composition(sequence):
    total = len(sequence)
    return {
        'A (%)': (sequence.count('A') / total) * 100,
        'T (%)': (sequence.count('T') / total) * 100,
        'G (%)': (sequence.count('G') / total) * 100,
        'C (%)': (sequence.count('C') / total) * 100,
    }

# Function to analyze sequences and store results in a Pandas DataFrame
def analyze_sequences(sequences):
    results = []

    for seq_id, sequence in sequences.items():
        seq_length = len(sequence)
        gc_content = calculate_gc_content(sequence)
        rev_complement = reverse_complement(sequence)
        transcription = transcribe_to_rna(sequence)
        composition = nucleotide_composition(sequence)

        # Append analysis to the results list
        results.append({
            'Sequence ID': seq_id,
            'Length': seq_length,
            'GC Content (%)': gc_content,
            'A (%)': round(composition['A (%)'], 2),
            'T (%)': round(composition['T (%)'], 2),
            'G (%)': round(composition['G (%)'], 2),
            'C (%)': round(composition['C (%)'], 2),
            'Reverse Complement (first 50)': rev_complement[:50],
            'RNA Transcription (first 50)': transcription[:50]
        })

    # Convert results to a DataFrame
    return pd.DataFrame(results)

# Function to calculate basic statistics using numpy
def calculate_statistics(df):
    gc_content = df['GC Content (%)'].to_numpy()
    lengths = df['Length'].to_numpy()

    stats = {
        'GC Content': {
            'Mean': np.mean(gc_content),
            'Median': np.median(gc_content),
            'Variance': np.var(gc_content),
            'Standard Deviation': np.std(gc_content),
        },
        'Lengths': {
            'Mean': np.mean(lengths),
            'Median': np.median(lengths),
            'Variance': np.var(lengths),
            'Standard Deviation': np.std(lengths),
        }
    }
    return stats

# Function to plot GC content distribution using Matplotlib
def plot_gc_content_distribution(df, output_path):
    plt.figure(figsize=(8, 6))
    plt.hist(df['GC Content (%)'], bins=10, color='blue', edgecolor='black', alpha=0.7)
    plt.title('GC Content Distribution')
    plt.xlabel('GC Content (%)')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.savefig(output_path)
    plt.close()

# Function to plot sequence length distribution using Matplotlib
def plot_sequence_length_distribution(df, output_path):
    plt.figure(figsize=(8, 6))
    plt.hist(df['Length'], bins=10, color='green', edgecolor='black', alpha=0.7)
    plt.title('Sequence Length Distribution')
    plt.xlabel('Length (bp)')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.savefig(output_path)
    plt.close()

# Function to generate a scatter plot for GC content vs. sequence length
def plot_gc_vs_length(df, output_path):
    plt.figure(figsize=(8, 6))
    plt.scatter(df['Length'], df['GC Content (%)'], c='purple', alpha=0.7, edgecolors='black')
    plt.title('GC Content vs. Sequence Length')
    plt.xlabel('Length (bp)')
    plt.ylabel('GC Content (%)')
    plt.grid(True)
    plt.savefig(output_path)
    plt.close()

# Function to generate a PDF report
def generate_pdf_report(df, stats, output_path, gc_plot, length_plot, scatter_plot):
    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()
    pdf.set_font("Arial", size=12)

    # Title
    pdf.set_font("Arial", size=16, style='B')
    pdf.cell(200, 10, txt="Sequence Analysis Report", ln=True, align='C')
    pdf.ln(10)

    # Summary statistics
    pdf.set_font("Arial", size=12)
    pdf.cell(200, 10, txt="--- Summary Statistics ---", ln=True)
    for key, value in stats.items():
        pdf.cell(200, 10, txt=f"{key}:", ln=True)
        for stat, result in value.items():
            pdf.cell(200, 10, txt=f"  {stat}: {result:.2f}", ln=True)
    pdf.ln(10)

    # Sequence details
    pdf.cell(200, 10, txt="--- Sequence Details ---", ln=True)
    for i, row in df.iterrows():
        pdf.cell(200, 10, txt=f"ID: {row['Sequence ID']}, Length: {row['Length']}, GC Content: {row['GC Content (%)']:.2f}%", ln=True)

    # Add plots
    pdf.add_page()
    pdf.cell(200, 10, txt="--- Plots ---", ln=True)
    pdf.ln(10)
    pdf.image(gc_plot, x=10, y=None, w=180)
    pdf.ln(85)
    pdf.image(length_plot, x=10, y=None, w=180)
    pdf.ln(85)
    pdf.image(scatter_plot, x=10, y=None, w=180)

    # Save PDF
    pdf.output(output_path)
    print(f"PDF report saved to {output_path}")

# Main function to run the program
def main():
    file_path = input("Enter the path to your FASTA file: ").strip()

    if not os.path.isfile(file_path):
        print("File not found. Please check the file path and try again.")
        return

    # Load sequences from the FASTA file
    sequences = load_fasta(file_path)
    if not sequences:
        print("No sequences found in the FASTA file.")
        return

    # Analyze sequences and store the results in a DataFrame
    df = analyze_sequences(sequences)

    # Calculate summary statistics
    stats = calculate_statistics(df)

    # Plot and save graphs
    gc_plot = "gc_content_distribution.png"
    length_plot = "sequence_length_distribution.png"
    scatter_plot = "gc_vs_length.png"

    plot_gc_content_distribution(df, gc_plot)
    plot_sequence_length_distribution(df, length_plot)
    plot_gc_vs_length(df, scatter_plot)

    # Generate PDF report
    output_pdf = "sequence_analysis_report.pdf"
    generate_pdf_report(df, stats, output_pdf, gc_plot, length_plot, scatter_plot)

if __name__ == "__main__":
    main()
