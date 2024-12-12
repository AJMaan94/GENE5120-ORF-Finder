import sequenceAnalysis
import argparse

class CommandLine:
    """
    Handles command line argument parsing for the ORF finder program.
    """
    def __init__(self, inOpts=None):
        self.parser = argparse.ArgumentParser(
            description='Program to analyze DNA sequences and find ORFs',
            epilog='Data will be distributed in the desired file.',
            add_help=True,
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
        )
        
        self.parser.add_argument('inFile', action='store', help='input file name')
        self.parser.add_argument('outFile', action='store', help='output file name')
        self.parser.add_argument(
            '-lG', '--longestGene',
            action='store',
            nargs='?',
            const=True,
            default=False,
            help='longest Gene in an ORF'
        )
        self.parser.add_argument(
            '-mG', '--minGene',
            type=int,
            choices=(0, 100, 200, 300, 500, 1000),
            default=100,
            action='store',
            help='minimum Gene length'
        )
        self.parser.add_argument(
            '-s', '--start',
            action='append',
            default=['ATG'],
            nargs='?',
            help='start Codon'
        )
        self.parser.add_argument(
            '-t', '--stop',
            action='append',
            default=['TAG', 'TGA', 'TAA'],
            nargs='?',
            help='stop Codon'
        )
        self.parser.add_argument(
            '-v', '--version',
            action='version',
            version='%(prog)s 0.1'
        )
        
        self.args = self.parser.parse_args(inOpts) if inOpts else self.parser.parse_args()


class Writer:
    """
    Handles writing ORF data to output files.
    
    Args:
        InfList: Dictionary containing ORF information
        outp: Name of output file
        h: Header information
    """
    def __init__(self, InfList, outp, h):
        self.InfList = InfList
        self.output = outp
        self.head = h
    
    def writeToFile(self):
        """
        Writes ORF data to the specified output file in a cleaned format.
        """
        with open(self.output, "a") as f:
            WriteDict = str(self.InfList)
            # Remove unnecessary characters
            chars_to_remove = "'{},:[]}{"
            for char in chars_to_remove:
                WriteDict = WriteDict.replace(char, "")
            
            writinglist = WriteDict.split("<")
            print("\n")
            f.write(self.head + "\n")
            
            for s in writinglist:
                f.write(s.lstrip() + "\n")
            f.write("\n")


def main(inCL=None):
    """
    Main function to process DNA sequences and find ORFs.
    
    Example usage:
    python findORFs.py inFile='tass2.fa' outFile='output.txt' --start=TTG --start=TGT 
    --stop=TAG --stop=TAA --minGene=100 --longestGene=True
    
    Note: Output data length is +2 compared to lab PDF as it counts the entire 
    stop codon as part of the reading frame.
    """
    # Initialize command line arguments
    myCommandLine = CommandLine(inCL)
    
    # Clean up input/output file names
    output = myCommandLine.args.outFile.replace("outFile=", "").replace("'", "")
    Input = myCommandLine.args.inFile.replace("inFile=", "").replace("'", "").replace(" ", "")
    
    # Initialize FastA reader
    myReader = sequenceAnalysis.FastAreader(Input)
    
    # Counter for unique nested dictionaries
    n = 0
    length = 0
    
    # Helper function for sorting by length
    def myL(e):
        return e[" Length "]
    
    # Process each sequence in the FASTA file
    for head, seq in myReader.readFasta():
        # Create ORF reader and process sequence
        TestFire = sequenceAnalysis.OrfReader(
            myCommandLine.args.start,
            myCommandLine.args.stop,
            head,
            seq,
            length,
            myCommandLine.args.minGene,
            n
        )
        InfoDict = TestFire.OrfAssembler()
        InfoList = list(InfoDict.values())
        InfoList.sort(reverse=True, key=myL)
        
        # Handle longest gene option
        if myCommandLine.args.longestGene in ["True", True]:
            InfoList = InfoList[:1]
        
        # Write results to file
        Write = Writer(InfoList, output, head)
        Write.writeToFile()
        
        # Update ORF counter
        n += len(InfoDict)
    
    print(f"Your data has been stored in {output}\n")
    print("Note! This code gives length as +2 compared to the output given on lab 5 "
          "as it chooses to count the entire stop codon as part of the reading frame. "
          "I believe this is more correct!\n")
    print("I added example user parameters in comments under main(inCl=None)")


if __name__ == "__main__":
    main()