import pysam
import argparse
from vcf_handler import VCFHandler

"""
files
/Users/jimin/CGL/nanopore_somatic/ONT_PMDV_vcf/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT_tumor_filtered.vcf.gz
/Users/jimin/CGL/nanopore_somatic/ONT_PMDV_vcf/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT_normal.vcf.gz
"""

parser = argparse.ArgumentParser()
parser.add_argument('--vcf', '-v', type=str, required=True, help='tumor vcf file path')
parser.add_argument('--vcftumorindex', '-i', type=str, help='tumor vcf index file path')
parser.add_argument('--vcfnormal', '-n', type=str, required=True, help='normal vcf file path')
parser.add_argument('--vcfnormalindex', '-z', type=str, help='normal vcf index file path')

parser.add_argument('--outfile', '-o', type=str, help='output file path')

# if contig argument not given, reads in whole vcf file
parser.add_argument('--contig', '-c', type=str, help='contig region to read in')

args = parser.parse_args()

def filter_out_normal(records, normal_pos_dict):
    for record in records:
        contig = record.chrom
        if record.pos not in normal_pos_dict[contig]:
            yield record

if __name__ == '__main__':
    tumor_vcf = args.vcf
    tumor_vcf_index = args.vcftumorindex
    tumor_vcf_h = VCFHandler(tumor_vcf, tumor_vcf_index)

    tumor_records = tumor_vcf_h.get_records()

    normal_vcf = args.vcfnormal
    normal_vcf_index = args.vcfnormalindex
    normal_vcf_h = VCFHandler(normal_vcf, normal_vcf_index)

    normal_records = normal_vcf_h.get_records()

    out_file_name = args.outfile
    out_file = pysam.VariantFile(out_file_name, 'w', header=tumor_vcf_h.header)

    tumor_pos_list = set()
    tumor_pos_list_filtered = set()
    normal_pos_list = set()
    normal_pos_dict = {}

    # get positions only from records
    chromosome = normal_vcf_h.get_record_attributes(normal_records, 'chrom')
    positions = normal_vcf_h.get_record_attributes(normal_records, 'pos')

    for c,p in zip(chromosome, positions):
        if c not in normal_pos_dict:
            normal_pos_dict[c] = set()
            normal_pos_dict[c].add(p)
        else:
            normal_pos_dict[c].add(p)


    filtered = filter_out_normal(tumor_records, normal_pos_dict)
    for record in filtered:
        out_file.write(record)




