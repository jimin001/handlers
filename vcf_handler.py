import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--vcf', '-v', type=str, required=True, help='vcf file path')
parser.add_argument('--vcfindex', '-i', type=str, required=True, help='vcf index file path')
args = parser.parse_args()


class SampleData:
    """
    Store information for GT:AP:GQ:DP:AD:VAF:REP.
    """
    def __init__(self, samples):
        self.samples = samples

        for sample in samples:
            self.sample = sample
            self.GT = self.sample[0]
            self.AP = self.sample[1]
            self.GQ = self.sample[2]
            self.DP = self.sample[3]
            self.AD = self.sample[4]
            self.VAF = self.sample[5]
            self.REP = self.sample[6]

class VCFHandler:
    """
    VCF Handler class.
    """
    def __init__(self, input_vcf_path, input_vcf_index_path):
        self.vcf = pysam.VariantFile(input_vcf_path, 'r', index_filename=input_vcf_index_path)
        self.header = self.vcf.header

    def get_header(self):
        """
        Get header of VCF.
        """
        return self.header

    def get_reads(self, contigs):
        """
        Read in vcf line by line in specified contig region
        """
        # returns iterator object
        return self.vcf.fetch(contig=contigs)

    def get_read_attributes(self, reads, attribute):
        """
        Get specified feature of vcf read.
        """
        # returns generator object of specified attribute
        for read in reads:
            if attribute == 'chrom':
                yield read.chrom
            elif attribute == 'pos':
                yield read.pos
            elif attribute == 'id':
                yield read.id
            elif attribute == 'ref':
                yield read.ref
            elif attribute == 'alt':
                yield read.alts
            elif attribute == 'qual':
                yield read.qual
            elif attribute == 'filter':
                yield read.filter.keys()
            elif attribute == 'info':
                yield read.info.keys()
            elif attribute == 'format':
                yield read.format.keys()
            elif attribute == 'sample':
                yield [y.values() for y in read.samples.itervalues()]

    def get_read_samples(self, samples_list, type):
        """
        Return generator object of specified data type of sample.
        :param samples_list: samples attribute generator object or list
        :param type: GT:AP:GQ:DP:AD:VAF:REP
        :return:
        """
        for s in samples_list:
            samples = SampleData(s)
            if type == 'GT':
                yield samples.GT
            elif type == 'AP':
                yield samples.AP
            elif type == 'GQ':
                yield samples.GQ
            elif type == 'DP':
                yield samples.DP
            elif type == 'AD':
                yield samples.AD
            elif type == 'VAF':
                yield samples.VAF
            elif type == 'REP':
                yield samples.REP
            else:
                return 'Incorrect sample type'

if __name__ == '__main__':
    vcf_file = args.vcf
    vcf_index = args.vcfindex

    vcf_handler = VCFHandler(vcf_file, vcf_index)
    vcf_header = vcf_handler.get_header()
    #print(vcf_header)

    reads = vcf_handler.get_reads(contigs='chr1')
    #for read in reads:
        #print(read)

    samples_list = vcf_handler.get_read_attributes(reads, 'sample')
    #print(samples_list)

    #print(vcf_handler.get_read_samples(samples_list, 'GT'))

    # all 'GT'
    #print([x for x in vcf_handler.get_read_samples(samples_list, 'GT')])











