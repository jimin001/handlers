import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--vcf', '-v', type=str, required=True, help='vcf file path')
parser.add_argument('--vcfindex', '-i', type=str, help='vcf index file path')
parser.add_argument('--outfile', '-o', type=str, help='output file path')
# if contig argument not given, reads in whole vcf file
parser.add_argument('--contig', '-c', type=str, help='contig region to read in')
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
    def __init__(self, input_vcf_path, input_vcf_index_path, contig):
        self.vcf = pysam.VariantFile(input_vcf_path, 'r', index_filename=input_vcf_index_path)
        self.header = self.vcf.header
        self.contig = contig

    def get_header(self):
        """
        Get header of VCF.
        """
        return self.header

    def get_records(self):
        """
        Read in vcf line by line in specified contig region
        """
        # returns iterator object
        return self.vcf.fetch(self.contig)

    def get_record_attributes(self, records, attribute):
        """
        Get specified feature of vcf record.
        """
        # returns generator object of specified attribute
        for record in records:
            if attribute == 'chrom':
                yield record.chrom
            elif attribute == 'pos':
                yield record.pos
            elif attribute == 'id':
                yield record.id
            elif attribute == 'ref':
                yield record.ref
            elif attribute == 'alt':
                yield record.alts
            elif attribute == 'qual':
                yield record.qual
            elif attribute == 'filter':
                yield record.filter.keys()
            elif attribute == 'info':
                yield record.info.keys()
            elif attribute == 'format':
                yield record.format.keys()
            elif attribute == 'sample':
                yield [y.values() for y in record.samples.itervalues()]

    def get_record_samples(self, samples_list, type):
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
    contig = args.contig

    vcf_handler = VCFHandler(vcf_file, vcf_index, contig)
    vcf_header = vcf_handler.get_header()
    #print(vcf_header)

    records = vcf_handler.get_records()
    for record in records:
        print(record)

    samples_list = vcf_handler.get_record_attributes(records, 'sample')
    #print(samples_list)

    #print(vcf_handler.get_record_samples(samples_list, 'GT'))

    # print all 'GT'
    #print([x for x in vcf_handler.get_record_samples(samples_list, 'GT')])











