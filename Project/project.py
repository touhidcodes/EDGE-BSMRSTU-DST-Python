# BioSeqAnalyzer: Biological Sequence Analysis and Visualization with Python

# import all Modules and Libraries
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from fpdf import FPDF

# Function to select FASTA files in the current directory
def selectFile():
    fastaFiles = [file for file in os.listdir() if file.endswith(".fasta")]

    if not fastaFiles:
        print("No FASTA files found in the current directory.")
        return
     
    # Display available FASTA files
    print("\nAvailable FASTA files:")
    for index, fileName in enumerate(fastaFiles, 1):
        print(f"{index}. {fileName}")
    
    userChoice = input("\nEnter the number of the FASTA file you want to use: ")

    # Validate user input
    if userChoice.isdigit():
        fileIndex = int(userChoice) - 1
        if 0 <= fileIndex and fileIndex < len(fastaFiles):
            selectedFile = fastaFiles[fileIndex]
            print(f"\nYou selected: {selectedFile}")
            return selectedFile
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

# Function to calculate GC content using BioPython
def calculateGCContent(sequence):
    return gc_fraction(sequence) * 100

# Function to get the reverse complement of a DNA sequence using BioPython
def reverseComplement(sequence):
    seqObj = Seq(sequence)
    return str(seqObj.reverse_complement())

# Function to transcribe DNA to RNA (T -> U)
def transcribeToRNA(sequence):
    seqObj = Seq(sequence)
    return str(seqObj.transcribe())

# Function to calculate nucleotide composition
def nucleotideComposition(sequence):
    total = len(sequence)
    return {
        'A (%)': (sequence.count('A') / total) * 100,
        'T (%)': (sequence.count('T') / total) * 100,
        'G (%)': (sequence.count('G') / total) * 100,
        'C (%)': (sequence.count('C') / total) * 100,
    }

# Function to analyze sequences and store results in a Pandas DataFrame
def analyzeSequences(sequences):
    results = []

    for seqId, sequence in sequences.items():
        results.append({
            'Sequence ID': seqId,
            'Length': len(sequence),
            'GC Content (%)': calculateGCContent(sequence),
            'A (%)': round(nucleotideComposition(sequence)['A (%)'], 2),
            'T (%)': round(nucleotideComposition(sequence)['T (%)'], 2),
            'G (%)': round(nucleotideComposition(sequence)['G (%)'], 2),
            'C (%)': round(nucleotideComposition(sequence)['C (%)'], 2),
            'Reverse Complement (first 50)': reverseComplement(sequence)[:50],
            'RNA Transcription (first 50)': transcribeToRNA(sequence)[:50]
        })
    return pd.DataFrame(results)

# Function to calculate basic statistics using numpy
def calculateStats(data):
    return {
        'Mean': np.mean(data),
        'Median': np.median(data),
        'Variance': np.var(data),
        'Standard Deviation': np.std(data)
    }

# Function to calculate GC  Content and Length Stats
def calculateBioStats(df):
    gcContent = df['GC Content (%)'].to_numpy()
    lengths = df['Length'].to_numpy()

    stats = {
        'GC Content': calculateStats(gcContent),
        'Lengths': calculateStats(lengths)
    }
    return stats

# Function to plot GC content distribution using Matplotlib
def plotGcContentDistribution(df, outputPath):
    plt.figure(figsize=(8, 6))
    plt.hist(df['GC Content (%)'], bins=15, color='skyblue', edgecolor='black', alpha=0.7)
    plt.title('GC Content Distribution', fontsize=16)
    plt.xlabel('GC Content (%)', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(outputPath)
    plt.close()

# Function to plot sequence length distribution using Matplotlib
def plotSequenceLengthDistribution(df, outputPath):
    plt.figure(figsize=(8, 6))
    plt.hist(df['Length'], bins=15, color='lightgreen', edgecolor='black', alpha=0.7)
    plt.title('Sequence Length Distribution', fontsize=16)
    plt.xlabel('Length (bp)', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(outputPath)
    plt.close()

# Function to generate a scatter plot for GC content vs. sequence length
def plotGcVsLength(df, outputPath):
    plt.figure(figsize=(8, 6))
    plt.scatter(df['Length'], df['GC Content (%)'], c='purple', alpha=0.7, edgecolors='black', s=100)
    plt.title('GC Content vs. Sequence Length', fontsize=16)
    plt.xlabel('Length (bp)', fontsize=14)
    plt.ylabel('GC Content (%)', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(outputPath)
    plt.close()

# Function to generate a Pair Plot for GC Content and Sequence Length
def plotPairPlot(df, outputPath):
    plt.figure(figsize=(8, 6))
    sns.pairplot(df[['GC Content (%)', 'Length']])
    plt.title('Pair Plot of GC Content and Length', fontsize=12)
    plt.tight_layout()
    plt.savefig(outputPath)
    plt.close()

# Function to generate a PDF report with plots
def generatePDFReport(stats, gcPlot, lengthPlot, scatterPlot, pairPlot, outputFolder):
    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()

    # Title of the PDF
    pdf.set_font('Arial', 'B', 16)
    pdf.cell(200, 10, txt=f"{outputFolder} Sequence Analysis Report", ln=True, align='C')
    pdf.ln(10) 

    # Add Summary Stats to the PDF
    pdf.set_font('Arial', '', 12)
    pdf.multi_cell(0, 10, txt=
                              f"GC Content Stats:\nMean: {stats['GC Content']['Mean']:.2f}\n"
                              f"Median: {stats['GC Content']['Median']:.2f}\n"
                              f"Variance: {stats['GC Content']['Variance']:.2f}\n"
                              f"Standard Deviation: {stats['GC Content']['Standard Deviation']:.2f}\n\n"
                              f"Sequence Length Stats:\nMean: {stats['Lengths']['Mean']:.2f}\n"
                              f"Median: {stats['Lengths']['Median']:.2f}\n"
                              f"Variance: {stats['Lengths']['Variance']:.2f}\n"
                              f"Standard Deviation: {stats['Lengths']['Standard Deviation']:.2f}\n")

    # Add the images of plots to the PDF
    pdf.add_page()
    pdf.cell(200, 10, txt="--- Analysis Plots ---", ln=True)
    pdf.ln(5)
    pdf.image(gcPlot, x=10, y=None, w=150)
    pdf.ln(10)
    pdf.image(lengthPlot, x=10, y=None, w=150)
    pdf.ln(10)
    pdf.image(scatterPlot, x=10, y=None, w=150)
    pdf.ln(10)
    pdf.image(pairPlot, x=10, y=None, w=120)

    # Save the PDF to a file
    pdfOutputPath = os.path.join(outputFolder, "analysis_report.pdf")
    pdf.output(pdfOutputPath)

    print(f"Analysis report generated and saved at: {pdfOutputPath}")

def main():
    
    #  Select FASTA file
    selectedFile = selectFile()
    
    # Load sequences from the FASTA file
    sequences = getSequences(selectedFile)
    if not sequences:
        print("No sequences found in the FASTA file.")
        return
    
    # Analyze sequences and store the results in a DataFrame
    df = analyzeSequences(sequences)

    # Calculate summary statistics
    stats = calculateBioStats(df)

    # Plot and save graphs to AnalysisResult folder
    outputFolder = os.path.splitext(selectedFile)[0]
    os.makedirs(outputFolder, exist_ok=True)

    # Set file paths for each plot
    gcPlot = os.path.join(outputFolder, "gc_content_distribution.png")
    lengthPlot = os.path.join(outputFolder, "sequence_length_distribution.png")
    scatterPlot = os.path.join(outputFolder, "gc_vs_length.png")
    pairPlot = os.path.join(outputFolder, "pair_plot_gc_length.png")

    # Generate and save plots
    plotGcContentDistribution(df, gcPlot)
    plotSequenceLengthDistribution(df, lengthPlot)
    plotGcVsLength(df, scatterPlot)
    plotPairPlot(df, pairPlot)

   # Generate PDF save analysis file
    generatePDFReport(stats, gcPlot, lengthPlot, scatterPlot, pairPlot, outputFolder)

# Run the program
if __name__ == "__main__":
    main()