import pysam
import argparse

# takes in input file
parser = argparse.ArgumentParser()
parser.add_argument('--bamFile', '-b', type=str, help='input file path of bam file')
parser.add_argument('--contig', '-c', type=str, help='contig region to read in')
parser.add_argument('--attribute', '-a', type=str, help='get specific attribute of reads')

args = parser.parse_args()

class BamHandler:
    def __init__(self, bamfile_path):
        self.bamfile_path = bamfile_path
        self.contig = None
        self.bamfile = pysam.AlignmentFile(bamfile_path, 'rb')
        self.header = self.bamfile.header

    def get_header(self):
        return self.header

    def get_reads(self):
        # return iterator object
        return self.bamfile.fetch(self.contig)

    def get_read_attributes(self, reads, attribute):
        for read in reads:
            if attribute == 'qname':
                yield read.query_name
            elif attribute == 'flag':
                yield read.flag
            elif attribute == 'rname':
                yield read.reference_id
            elif attribute == 'pos':
                yield read.reference_start
            elif attribute == 'mapq':
                yield read.mapping_quality
            elif attribute == 'cigar':
                yield read.cigarstring
            elif attribute == 'rnext':
                yield read.next_reference_id
            elif attribute == 'pnext':
                yield  read.next_reference_start
            elif attribute == 'tlen':
                yield read.template_length
            elif attribute == 'seq':
                yield read.query_sequence
            elif attribute == 'qual':
                # array type
                yield read.query_qualities
            elif attribute == 'query length':
                yield read.query_length
            else:
                return 'Unknown attribute entered'


if __name__ == '__main__':
    inFile = args.bamFile
    contig = args.contig
    attribute = args.attribute

    bam_handler = BamHandler(inFile)

    if contig:
        bam_handler.contig = contig

    reads = bam_handler.get_reads()

    if not attribute:
        for read in reads:
            print(read)

    if attribute:
        attribute_reads = bam_handler.get_read_attributes(reads, attribute)

        for read in attribute_reads:
            print(read)

    #position = bam_handler.get_read_attributes(reads, 'qual')
    #print(position)
    #print(bam_handler.header)




