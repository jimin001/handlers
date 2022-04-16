import pysam
import argparse
from vcf_handler import VCFHandler

parser = argparse.ArgumentParser()
parser.add_argument('--vcf', '-v', type=str, required=True, help='vcf file path')
parser.add_argument('--vcfindex', '-i', type=str, required=True, help='vcf index file path')
parser.add_argument('--outfile', '-o', type=str, help='output file path')
parser.add_argument('--filter', '-f', type=str, help='filter type to filter out')

args = parser.parse_args()

def filter_by_filter(records, filter_type):
    """
    Filter out filter types given by user.
    """
    for record in records:
        if filter_type not in record.filter.keys():
            yield record


if __name__ == '__main__':
    vcf_file = args.vcf
    vcf_index = args.vcfindex
    vcf_h = VCFHandler(vcf_file, vcf_index)

    out_file_name = args.outfile
    out_file = pysam.VariantFile(out_file_name, 'w', header=vcf_h.header)

    filter_type = args.filter

    records = vcf_h.get_records()

    filtered = filter_by_filter(records, filter_type)

    # write to vcf
    for record in filtered:
        out_file.write(record)



