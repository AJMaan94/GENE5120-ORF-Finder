#!/usr/bin/env python3
import sys

class FastAreader:
    """Handles reading of FASTA format files."""
    
    def __init__(self, fname=''):
        """Constructor: saves attribute fname"""
        self.fname = fname
    
    def doOpen(self):
        """Opens the file or returns stdin if no filename provided."""
        return sys.stdin if self.fname == '' else open(self.fname)
    
    def readFasta(self, SeqStart='>'):
        """
        Reads a FASTA file and yields tuples of header and sequence.
        
        Args:
            SeqStart (str): Character that indicates the start of a new sequence (default: '>')
        
        Yields:
            tuple: (header, sequence) pairs from the FASTA file
        """
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
            
            # Skip to first FASTA header
            line = fileH.readline()
            while not line.startswith(SeqStart):
                line = fileH.readline()
            header = line[1:].rstrip()
            
            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()
            
            yield header, sequence


class OrfReader:
    """
    Handles reading and analyzing Open Reading Frames (ORFs) in DNA sequences.
    
    Args:
        SCSet: Set of start codons
        ECSet: Set of end codons
        Header: Sequence header
        Sequence: DNA sequence
        Length: Gene length
        mini: Minimum length requirement
        number: ORF number
    """
    
    def __init__(self, SCSet, ECSet, Header, Sequence, Length, mini, number):
        # Initialize codon sets
        self.n = number
        self.SCSet = SCSet
        self.ECSet = ECSet
        self.geneLength = Length
        self.mini = mini
        
        # Create reverse codon sets
        self.revSetS = set()
        self.revSetE = set()
        self.RevSCSet = self.SReverseCodon(self.SCSet)
        self.RevECSet = self.EReverseCodon(self.ECSet)
        
        # Store sequence information
        self.Header = Header
        self.Sequence = Sequence
        
        # Initialize dictionaries
        self.FOrfDic = {}
        self.ROrfDic = {}
        self.OrfDic = {}
        self.TempFOrfDic = {}
        self.TempROrfDic = {}
    
    def OrfAssembler(self):
        """
        Assembles ORFs from both forward and reverse strands.
        
        Returns:
            dict: Combined dictionary of forward and reverse ORFs
        """
        for x in range(0, 3):
            self.FOrfDic = self.FindFSCodons(x)
            self.ROrfDic = self.FindRSCodons(x)
            if self.FOrfDic:
                self.OrfDic.update(self.FOrfDic)
            if self.ROrfDic:
                self.OrfDic.update(self.ROrfDic)
        return self.OrfDic
    
    def SReverseCodon(self, CodonSetS):
        """Creates reverse complement of start codons."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        for codon in CodonSetS:
            rev_codon = ''.join(complement[base] for base in reversed(codon))
            self.revSetS.add(rev_codon)
        return self.revSetS
    
    def EReverseCodon(self, CodonSetE):
        """Creates reverse complement of end codons."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        for codon in CodonSetE:
            rev_codon = ''.join(complement[base] for base in reversed(codon))
            self.revSetE.add(rev_codon)
        return self.revSetE
    
    def FindFSCodons(self, X):
        """
        Finds ORFs in forward direction for a given frame.
        
        Args:
            X: Reading frame offset (0, 1, or 2)
        
        Returns:
            dict: Dictionary of ORFs found in this frame
        """
        self.framef = X
        for i in range(self.framef, len(self.Sequence) + 1, 3):
            codon = self.Sequence[i:i+3]
            if codon in self.SCSet:
                self.StartPos = i + self.geneLength + 1
                self.EndPos = self.endCFinder(i) + self.geneLength + 3
                self.ReadingFrame = self.framef + 1
                self.Length = self.EndPos - self.StartPos + 3
                
                if self.Length >= self.mini:
                    self.n += 1
                    self.TempFOrfDic[self.n] = {
                        "Start Position: ": self.StartPos,
                        " End Position": self.EndPos,
                        " Reading Frame ": self.ReadingFrame,
                        " Length ": self.Length,
                        "<": ""
                    }
                    return self.TempFOrfDic
        return {}
    
    def FindRSCodons(self, X):
        """
        Finds ORFs in reverse direction for a given frame.
        
        Args:
            X: Reading frame offset (0, 1, or 2)
        
        Returns:
            dict: Dictionary of ORFs found in this frame
        """
        self.frameR = X
        for y in range(len(self.Sequence) - 3 - self.frameR, -1, -3):
            codon = self.Sequence[y:y+3]
            if codon in self.RevSCSet:
                self.REndPos = self.revendCFinder(y) + 1 + self.geneLength + 3
                self.RReadingFrame = 0 - self.frameR - 1
                self.RStartPos = y + 1 + self.geneLength
                self.Length = self.RStartPos - self.REndPos + 3
                
                if self.Length >= self.mini:
                    self.n += 1
                    self.TempROrfDic[self.n] = {
                        "Start Position ": self.RStartPos,
                        " End Position": self.REndPos,
                        " Reading Frame ": self.RReadingFrame,
                        " Length ": self.Length,
                        "<": ""
                    }
                    return self.TempROrfDic
        return {}
    
    def endCFinder(self, startPosition):
        """Finds the next stop codon in forward direction."""
        for o in range(startPosition, len(self.Sequence) + 1, 3):
            codon = self.Sequence[o:o+3]
            if codon in self.ECSet and o <= len(self.Sequence):
                return o
        return len(self.Sequence)
    
    def revendCFinder(self, startPosition):
        """Finds the next stop codon in reverse direction."""
        for p in range(startPosition, -1, -3):
            codon = self.Sequence[p:p+3]
            if codon in self.RevECSet and p >= 0:
                return p
        return 0


class ProteinParam:
    """
    Calculates various physical and chemical properties of a protein sequence.
    """
    
    # Amino acid property tables
    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }
    mwH2O = 18.015
    aa2abs280 = {'Y': 1490, 'W': 5500, 'C': 125, 'F': 200}
    aa2chargePos = {'K': 10.5, 'R': 12.4, 'H': 6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34
    
    def __init__(self, protein):
        """Initialize with a protein sequence."""
        self.AminoList = list(self.aa2mw.keys())
        self.Chain = protein.upper()
        self.StringLength = len(protein)
        self.AminoCount = None
        self.AADict = None
        self.aaCount()  # Initialize amino acid counts
    
    def aaCount(self):
        """Count the number of each amino acid in the sequence."""
        self.AminoCount = [self.Chain.count(aa) for aa in self.AminoList]
        self.AADict = dict(zip(self.AminoList, self.AminoCount))
        return sum(self.AminoCount)
    
    def pI(self):
        """Calculate the isoelectric point."""
        return self._charge_()
    
    def aaComposition(self):
        """Return the amino acid composition dictionary."""
        return self.AADict
    
    def _charge_(self):
        """Calculate the net charge at different pH values."""
        best_ph = 0
        lowest_charge_diff = float('inf')
        
        for n in range(1400):
            ph = n * 0.01
            pos_charge = 0.0
            neg_charge = 0.0
            
            # Calculate positive charges
            for aa, pka in self.aa2chargePos.items():
                count = self.AADict.get(aa, 0)
                pos_charge += count * (10**pka) / (10**pka + 10**ph)
            pos_charge += (10**self.aaNterm) / (10**self.aaNterm + 10**ph)
            
            # Calculate negative charges
            for aa, pka in self.aa2chargeNeg.items():
                count = self.AADict.get(aa, 0)
                neg_charge += count * (10**ph) / (10**pka + 10**ph)
            neg_charge += (10**ph) / (10**self.aaCterm + 10**ph)
            
            charge_diff = abs(pos_charge - neg_charge)
            if charge_diff < lowest_charge_diff:
                lowest_charge_diff = charge_diff
                best_ph = ph
        
        return best_ph
    
    def molarExtinction(self):
        """Calculate the molar extinction coefficient."""
        return sum(count * self.aa2abs280[aa] 
                  for aa, count in self.AADict.items() 
                  if aa in self.aa2abs280)
    
    def massExtinction(self):
        """Calculate the mass extinction coefficient."""
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0
    
    def molecularWeight(self):
        """Calculate the molecular weight of the protein."""
        return (self.mwH2O + 
                sum(count * (self.aa2mw[aa] - self.mwH2O) 
                    for aa, count in self.AADict.items()))


class NucParams:
    """
    Analyzes nucleotide sequences for various parameters including codon usage,
    amino acid composition, and nucleotide composition.
    """
    
    def __init__(self):
        self.Chain = ""
        self.CodonCompD = {key: 0 for key in rnaCodonTable}
        self.aaCompD = {value: 0 for value in set(rnaCodonTable.values())}
        self.NucD = {"A": 0, "G": 0, "C": 0, "T": 0, "U": 0, "N": 0}
        self.AllParameters = {}
    
    def addSequence(self, filename):
        """Process sequences from a FASTA file."""
        myReader = FastAreader(filename)
        
        for head, seq in myReader.readFasta():
            self.ChainUpper = seq.upper()
            self.nucTable = self.nucComposition()
            self.codonTable = self.codonComposition()
        
        self.aaTable = self.aaComposition()
        nuc_count = self.nucCount()
        
        self.AllParameters = {
            "codonComposition": self.codonTable,
            "aaComposition": self.aaTable,
            "nucComposition": self.nucTable,
            "Nuc Count": nuc_count
        }
        
        return self.AllParameters
    
    def aaComposition(self):
        """Calculate amino acid composition."""
        for codon, aa in rnaCodonTable.items():
            self.aaCompD[aa] = self.aaCompD[aa] + self.CodonCompD[codon]
        return self.aaCompD
    
    def nucComposition(self):
        """Calculate nucleotide composition."""
        for key in self.NucD:
            self.NucD[key] = self.NucD[key] + self.ChainUpper.count(key)
        return self.NucD
    
    def codonComposition(self):
        """Calculate codon composition."""
        for i in range(0, len(self.ChainUpper), 3):
            codon = self.ChainUpper[i:i+3].replace("T", "U")
            if "N" not in codon and codon in rnaCodonTable:
                self.CodonCompD[codon] += 1
        return self.CodonCompD
    
    def nucCount(self):
        """Return total nucleotide count."""
        return sum(self.NucD.values())