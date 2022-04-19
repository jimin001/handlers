import pysam
import argparse

# takes in input files
parser = argparse.ArgumentParser()
parser.add_argument('--fasta', '-f', type=str, required=True, help='input fasta file path')
parser.add_argument('--index', '-i', type=str, help='input fasta index file path')
parser.add_argument('--references', '-x', action='store_true', help='will return list of references if flag given')
parser.add_argument('--region', '-r', type=str, help='[optional] specify region of fasta file')

args = parser.parse_args()

"""
python3 fasta_handler.py 
-f /Users/jimin/PEPPER_clr/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta 
-r chr1:1-1100 
-x
"""

class FastaHandler:
    def __init__(self, input_fasta_path):
        self.input_index_path = None
        self.region = None
        self.fasta = pysam.FastaFile(input_fasta_path, self.input_index_path)

    def get_sequences(self):
        return self.fasta.fetch(region=self.region)

    def get_references(self):
        return self.fasta.references


if __name__ == '__main__':
    inFile = args.fasta
    indexFile = args.index
    region = args.region
    ref = args.references

    fasta_handler = FastaHandler(inFile)

    if indexFile:
        fasta_handler.input_index_path = indexFile

    if region:
        fasta_handler.region = region

    if ref:
        references = fasta_handler.get_references()
        print(references)

    sequences = fasta_handler.get_sequences()
    print(sequences)





