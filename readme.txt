ORF Finder and Sequence Analysis Tool
Version 0.1
Overview
This tool provides functionality for finding Open Reading Frames (ORFs) in DNA sequences and analyzing various properties of nucleotide and protein sequences. It consists of two main Python scripts:

findORFs.py: Main script for ORF finding
sequenceAnalysis.py: Module containing sequence analysis functionality

Requirements

Python 3.6 or higher
Input files in FASTA format

Installation

Place both Python files in the same directory:
Copyyour_directory/
├── findORFs.py
└── sequenceAnalysis.py

Make sure the scripts have executable permissions:
bashCopychmod +x findORFs.py


Usage
Basic Usage
bashCopypython findORFs.py inFile='input.fa' outFile='output.txt'
Command Line Arguments

inFile: Input FASTA file (required)
outFile: Output file name (required)
-lG, --longestGene: Only output the longest gene in each ORF (optional)
-mG, --minGene: Minimum gene length (choices: 0,100,200,300,500,1000; default: 100)
-s, --start: Start codon(s) (default: ['ATG'])
-t, --stop: Stop codon(s) (default: ['TAG','TGA','TAA'])
-v, --version: Show program version

Extended Usage Examples

Basic search with default parameters:

bashCopypython findORFs.py inFile='sequence.fa' outFile='results.txt'

Search with custom minimum gene length:

bashCopypython findORFs.py inFile='sequence.fa' outFile='results.txt' --minGene=300

Search with additional start codons:

bashCopypython findORFs.py inFile='sequence.fa' outFile='results.txt' --start=TTG --start=GTG

Only output longest genes:

bashCopypython findORFs.py inFile='sequence.fa' outFile='results.txt' --longestGene=True

Full example with multiple parameters:

bashCopypython findORFs.py inFile='sequence.fa' outFile='results.txt' --start=TTG --start=GTG --stop=TAG --stop=TAA --minGene=100 --longestGene=True
Output Format
The output file contains the following information for each ORF:

Start Position
End Position
Reading Frame
Length

Example output:
Copy>Sequence_Header
Start Position: 1234 End Position: 2345 Reading Frame: 1 Length: 1112
Start Position: 3456 End Position: 4567 Reading Frame: -2 Length: 1112
Notes

The program counts the entire stop codon as part of the reading frame, resulting in lengths that are +2 compared to some other tools
Reading frames are numbered 1, 2, 3 for forward strand and -1, -2, -3 for reverse strand
All input sequences are automatically converted to uppercase
Invalid codons or sequences containing 'N' are skipped

Error Handling

Invalid FASTA files will generate an error message
Invalid command line arguments will display usage information
The program will create the output file if it doesn't exist, or append to it if it does

Additional Features
The sequenceAnalysis.py module provides additional functionality:

Protein parameter calculations (molecular weight, pI, etc.)
Nucleotide composition analysis
Codon usage analysis
Amino acid composition analysis

Limitations

Input files must be in proper FASTA format
Very large sequences may require significant processing time
Memory usage scales with sequence length

Support
For issues or questions, please create an issue on the repository or contact the developer:

Author: Amarjot Maan
Version: 0.1