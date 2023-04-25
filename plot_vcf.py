# plot different features of a VCF file

import pysam
import argparse
import matplotlib.pyplot as plt
import csv
import matplotlib.patches as mplpatches

"""
Running Command:
python3 plot_vcf.py \
--queryvcf query.vcf.gz \
--queryvcfindex query.vcf.gz.tbi \
--truthvcf truth.vcf.gz \
--truthvcfindex truth.vcf.gz.tbi \
--features sompy.features.csv \
--querytype 'clairs' \
--plot 'chromosome' \
--samplename 'name of sample' \
-o output_file_name.png
"""

parser = argparse.ArgumentParser()
parser.add_argument('--queryvcf', '-q', type=str, required=True, help='vcf file path')
parser.add_argument('--queryvcfindex', '-i', type=str, help='vcf index file path')
parser.add_argument('--truthvcf', '-t', type=str, required=True, help='vcf file path')
parser.add_argument('--truthvcfindex', '-x', type=str, help='vcf index file path')
parser.add_argument('--querytype', '-z', type=str, required=True, help='specify which variant caller query VCF is from, options: clairs, deepsomatic')
parser.add_argument('--features', '-f', type=str, required=True, help='features.csv from som.py')
parser.add_argument('--outfile', '-o', type=str, help='output file path')
parser.add_argument('--plot', '-p', type=str, required=True, help='specify which plot type, options: chromosome, histogram')
parser.add_argument('--samplename', '-s', type=str, help='sample name to label plot')
parser.add_argument('--minvaf', '-m', type=float, help='min vaf')
args = parser.parse_args()

class ClairSampleData:
    """
    Store information for GT:GQ:DP:AF:NAF:NDP:AU:CU:GU:TU.
    """
    def __init__(self, samples):
        self.samples = samples

        for sample in samples:
            self.sample = sample
            self.GT = self.sample[0]
            self.GQ = self.sample[1]
            self.DP = self.sample[2]
            self.AF = self.sample[3]
            self.NAF = self.sample[4]
            self.NDP = self.sample[5]
            self.AU = self.sample[6]
            self.CU = self.sample[7]
            self.GU = self.sample[8]
            self.TU = self.sample[9]

class DeepSomaticSampleData:
    """
    Store information for GT:AP:GQ:DP:AD:VAF:REP.
    """
    def __init__(self, samples):
        self.samples = samples

        for sample in samples:
            self.sample = sample
            self.GT = self.sample[0]
            self.GQ = self.sample[1]
            self.DP = self.sample[2]
            self.AD = self.sample[3]
            self.VAF = self.sample[4]
            self.PL = self.sample[5]

class TruthSampleData:
    """
    Store information for
    """
    def __init__(self, dict_info):
        self.info = dict_info

        for sample in dict_info:
            self.sample = sample
            if 'NVAF' in self.info:
                self.NVAF = self.info['NVAF']
            else:
                self.NVAF = None

            if 'TVAF' in self.info:
                self.TVAF = self.info['TVAF']
            else:
                self.TVAF = None

class VCFDictionary:
    """
    Read in CSV for som.py features file
    """
    def __init__(self, query_vcf, query_index, truth_vcf, truth_index):
        self.features_dict = {}
        self.qvcf = pysam.VariantFile(query_vcf, 'r', index_filename=query_index)
        self.qheader = self.qvcf.header
        self.qcontig = None
        self.tvcf = pysam.VariantFile(truth_vcf, 'r', index_filename=truth_index)
        self.theader = self.tvcf.header
        self.tcontig = None

    def sompyFeaturesReader(self, file_name):
        with open(file_name) as csvfile:
            csv_reader = csv.reader(csvfile)
            for row in csv_reader:
                chrom = row[1]
                pos = row[2]
                chrom_pos = str(chrom) + '_' + str(pos)
                classification = row[3]
                if pos != 'POS':
                    self.features_dict[chrom_pos] = [classification]

    def queryReader(self, records, querytype):
        for record in records:
            sample_info = [y.values() for y in record.samples.itervalues()]
            if querytype == 'clairs':
                sample_obj = ClairSampleData(sample_info)
                vaf = sample_obj.AF
            elif querytype == 'deepsomatic':
                sample_obj = DeepSomaticSampleData(sample_info)
                vaf = sample_obj.VAF[0]
                dp = sample_obj.DP
                ad = sample_obj.AD

            chrom = str(record.chrom)
            pos = str(record.pos)
            chrom_pos = chrom + '_' + pos

            if chrom == 'chrX':
                continue
            elif chrom == 'chrY':
                continue
            elif chrom == 'chrM':
                continue
            elif chrom_pos in self.features_dict:
                self.features_dict[chrom_pos] += [vaf, "query"]
                self.features_dict[chrom_pos] += [dp, ad]
            else:
                self.features_dict[chrom_pos] = [vaf, "query"]
                self.features_dict[chrom_pos] += [dp, ad]

    def truthReader(self, records):
        for record in records:
            sample_obj = TruthSampleData(record.info)
            chrom = str(record.chrom)
            pos = str(record.pos)
            chrom_pos = chrom + '_' + pos
            vaf = sample_obj.TVAF
            if vaf is not None:
                if chrom_pos in self.features_dict:
                    self.features_dict[chrom_pos] += [vaf, "truth"]
                else:
                    self.features_dict[chrom_pos] = [vaf, "truth"]


    def get_truth_header(self):
        """
        Get header of VCF.
        """
        return self.theader

    def get_query_header(self):
        """
        Get header of VCF.
        """
        return self.qheader

    def get_truth_records(self):
        """
        Read in vcf line by line in specified contig region
        """
        # returns iterator object
        return self.tvcf.fetch(self.tcontig)

    def get_query_records(self):
        """
        Read in vcf line by line in specified contig region
        """
        # returns iterator object
        return self.qvcf.fetch(self.qcontig)

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



def plotByChrom(vcf_dictionary, chrom, panel):
    count_FN = 0
    count_FP = 0
    count_TP = 0
    for position in vcf_dictionary:
        if len(vcf_dictionary[position]) >= 3:
            #print(position, vcf_dictionary[position])
            classification = vcf_dictionary[position][0]
            dict_chrom = position.split('_')[0]
            pos = int(position.split('_')[1])
            #print(classification, dict_chrom, pos)
            if dict_chrom == chrom:
                if classification == 'TP':
                    truth_vaf = float(vcf_dictionary[position][1])
                    query_vaf = float(vcf_dictionary[position][3])
                    # print(position, vcf_dictionary[position])
                    count_TP += 1
                    panel.plot(pos, truth_vaf,
                                    marker='s',
                                    markersize=2,
                                    linewidth=0,
                                    markeredgewidth=0,
                                    color='gray'
                                    )
                    panel.plot(pos, query_vaf,
                                    marker='s',
                                    markersize=2,
                                    linewidth=0,
                                    markeredgewidth=0,
                                    color='black'
                                    )

                elif classification == 'FN':
                    vaf = vcf_dictionary[position][1]
                    count_FN += 1
                    panel.plot(pos, vaf,
                                    marker='o',
                                    markersize=2,
                                    linewidth=0,
                                    markeredgewidth=0,
                                    color='blue'
                                    )
                elif classification == 'FP':
                    vaf = vcf_dictionary[position][1]
                    count_FP += 1
                    panel.plot(pos, vaf,
                                    marker='o',
                                    markersize=2,
                                    linewidth=0,
                                    markeredgewidth=0,
                                    color='red'
                                    )
    #print(count_FN, count_FP, count_TP)
def plotHistograms(vcf_dictionary, TP_truth_panel, TP_query_panel, FN_panel, FP_panel):
    TP_truth = []
    TP_query = []
    FN = []
    FP = []
    for position in vcf_dictionary:
        if len(vcf_dictionary[position]) >= 3:
            #print(position, vcf_dictionary[position])
            classification = vcf_dictionary[position][0]
            if classification == 'TP':
                truth_vaf = float(vcf_dictionary[position][1])
                query_vaf = float(vcf_dictionary[position][3])
                TP_truth.append(truth_vaf)
                TP_query.append(query_vaf)
            elif classification == 'FN':
                vaf = vcf_dictionary[position][1]
                FN.append(vaf)
            elif classification == 'FP':
                vaf = vcf_dictionary[position][1]
                FP.append(vaf)
    TP_truth_panel.hist(TP_truth, bins=50, color='gray')
    TP_query_panel.hist(TP_query, bins=50, color='black')
    FN_panel.hist(FN, bins=50, color='blue')
    FP_panel.hist(FP, bins=50, color='red')

def noTicks(panel):
    panel.tick_params(bottom=True,
                labelbottom=True,
                left=False,
                labelleft=False,
                right=False,
                labelright=False,
                top=False,
                labeltop=False)

if __name__ == '__main__':
    q_vcf_file = args.queryvcf
    q_vcf_index = args.queryvcfindex
    t_vcf_file = args.truthvcf
    t_vcf_index = args.truthvcfindex
    features_csv = args.features
    out_file_name = args.outfile
    query_type = args.querytype
    plot_type = args.plot
    sample_name = args.samplename

    # create VCFDictionary object
    vcf_object = VCFDictionary(q_vcf_file, q_vcf_index, t_vcf_file, t_vcf_index)

    # get records
    truth_records = vcf_object.get_truth_records()
    query_records = vcf_object.get_query_records()

    # read in data
    vcf_object.sompyFeaturesReader(features_csv)
    vcf_object.truthReader(truth_records)
    vcf_object.queryReader(query_records, query_type)

    # dictionary object that holds organized VCF information
    vcf_dictionary = vcf_object.features_dict


    # figure details
    figureWidth = 18
    figureHeight = 14

    plt.figure(figsize=(figureWidth, figureHeight))

    def chromPanels():
        # main panel dimensions
        chromWidth = 3
        chromHeight = 2

        legendWidth = 2
        legendHeight = 2

        panelWidth = 15
        panelHeight = 13

        # plot panels
        chromPanel21 = plt.axes([0.05, 0.05, chromWidth / figureWidth, chromHeight / figureHeight])
        chromPanel22 = plt.axes([0.24, 0.05, chromWidth / figureWidth, chromHeight / figureHeight])
        legendPanel = plt.axes([0.43, 0.05, legendWidth / figureWidth, legendHeight / figureHeight])

        chromPanel16 = plt.axes([0.05, 0.25, chromWidth / figureWidth, chromHeight / figureHeight])
        chromPanel17 = plt.axes([0.24, 0.25, chromWidth / figureWidth, chromHeight / figureHeight])
        chromPanel18 = plt.axes([0.43, 0.25, chromWidth / figureWidth, chromHeight / figureHeight])
        chromPanel19 = plt.axes([0.62, 0.25, chromWidth / figureWidth, chromHeight / figureHeight])
        chromPanel20 = plt.axes([0.81, 0.25, chromWidth / figureWidth, chromHeight / figureHeight])

        chromPanel11 = plt.axes([0.05, 0.45, chromWidth / figureWidth, chromHeight / figureHeight])
        chromPanel12 = plt.axes([0.24, 0.45, chromWidth / figureWidth, chromHeight / figureHeight])
        chromPanel13 = plt.axes([0.43, 0.45, chromWidth / figureWidth, chromHeight / figureHeight])
        chromPanel14 = plt.axes([0.62, 0.45, chromWidth / figureWidth, chromHeight / figureHeight])
        chromPanel15 = plt.axes([0.81, 0.45, chromWidth / figureWidth, chromHeight / figureHeight])

        chromPanel6 = plt.axes([0.05, 0.65, chromWidth / figureWidth, chromHeight / figureHeight])
        chromPanel7 = plt.axes([0.24, 0.65, chromWidth / figureWidth, chromHeight / figureHeight])
        chromPanel8 = plt.axes([0.43, 0.65, chromWidth / figureWidth, chromHeight / figureHeight])
        chromPanel9 = plt.axes([0.62, 0.65, chromWidth / figureWidth, chromHeight / figureHeight])
        chromPanel10 = plt.axes([0.81, 0.65, chromWidth / figureWidth, chromHeight / figureHeight])

        chromPanel1 = plt.axes([0.05, 0.85, chromWidth / figureWidth, chromHeight / figureHeight])
        chromPanel2 = plt.axes([0.24, 0.85, chromWidth / figureWidth, chromHeight / figureHeight])
        chromPanel3 = plt.axes([0.43, 0.85, chromWidth / figureWidth, chromHeight / figureHeight])
        chromPanel4 = plt.axes([0.62, 0.85, chromWidth / figureWidth, chromHeight / figureHeight])
        chromPanel5 = plt.axes([0.81, 0.85, chromWidth / figureWidth, chromHeight / figureHeight])

        chromPanel1.set_xlabel('chr1', fontsize=16)
        chromPanel2.set_xlabel('chr2', fontsize=16)
        chromPanel3.set_xlabel('chr3', fontsize=16)
        chromPanel4.set_xlabel('chr4', fontsize=16)
        chromPanel5.set_xlabel('chr5', fontsize=16)
        chromPanel6.set_xlabel('chr6', fontsize=16)
        chromPanel7.set_xlabel('chr7', fontsize=16)
        chromPanel8.set_xlabel('chr8', fontsize=16)
        chromPanel9.set_xlabel('chr9', fontsize=16)
        chromPanel10.set_xlabel('chr10', fontsize=16)
        chromPanel11.set_xlabel('chr11', fontsize=16)
        chromPanel12.set_xlabel('chr12', fontsize=16)
        chromPanel13.set_xlabel('chr13', fontsize=16)
        chromPanel14.set_xlabel('chr14', fontsize=16)
        chromPanel15.set_xlabel('chr15', fontsize=16)
        chromPanel16.set_xlabel('chr16', fontsize=16)
        chromPanel17.set_xlabel('chr17', fontsize=16)
        chromPanel18.set_xlabel('chr18', fontsize=16)
        chromPanel19.set_xlabel('chr19', fontsize=16)
        chromPanel20.set_xlabel('chr20', fontsize=16)
        chromPanel21.set_xlabel('chr21', fontsize=16)
        chromPanel22.set_xlabel('chr22', fontsize=16)

        chromPanel1.set_ylabel('VAF')
        chromPanel6.set_ylabel('VAF')
        chromPanel11.set_ylabel('VAF')
        chromPanel16.set_ylabel('VAF')
        chromPanel21.set_ylabel('VAF')

        circle1 = mplpatches.Circle((0.15, 0.8), radius=.05,
                                          facecolor='blue',
                                          edgecolor='black',
                                          linewidth=0)
        circle2 = mplpatches.Circle((0.15, 0.6), radius=.05,
                                          facecolor='red',
                                          edgecolor='black',
                                          linewidth=0)
        rectangle3 = mplpatches.Rectangle((0.1, 0.35), .1, .1,
                                          facecolor='gray',
                                          edgecolor='black',
                                          linewidth=0)
        rectangle4 = mplpatches.Rectangle((0.1, 0.15), .1, .1,
                                          facecolor='black',
                                          edgecolor='black',
                                          linewidth=0)
        legendPanel.add_patch(circle1)
        legendPanel.add_patch(circle2)
        legendPanel.add_patch(rectangle3)
        legendPanel.add_patch(rectangle4)

        legendPanel.text(0.26, 0.75, 'FN', fontsize=16, va='bottom', ha='left')
        legendPanel.text(0.26, 0.55, 'FP', fontsize=16, va='bottom', ha='left')
        legendPanel.text(0.26, 0.35, 'TP truth', fontsize=16, va='bottom', ha='left')
        legendPanel.text(0.26, 0.15, 'TP query', fontsize=16, va='bottom', ha='left')

        chromPanel1.tick_params(bottom=True,
                              labelbottom=True,
                              left=True,
                              labelleft=True,
                              right=False,
                              labelright=False,
                              top=False,
                              labeltop=False)
        chromPanel6.tick_params(bottom=True,
                                labelbottom=True,
                                left=True,
                                labelleft=True,
                                right=False,
                                labelright=False,
                                top=False,
                                labeltop=False)
        chromPanel11.tick_params(bottom=True,
                                labelbottom=True,
                                left=True,
                                labelleft=True,
                                right=False,
                                labelright=False,
                                top=False,
                                labeltop=False)
        chromPanel16.tick_params(bottom=True,
                                labelbottom=True,
                                left=True,
                                labelleft=True,
                                right=False,
                                labelright=False,
                                top=False,
                                labeltop=False)
        chromPanel21.tick_params(bottom=True,
                                labelbottom=True,
                                left=True,
                                labelleft=True,
                                right=False,
                                labelright=False,
                                top=False,
                                labeltop=False)
        legendPanel.tick_params(bottom=False,
                                 labelbottom=False,
                                 left=False,
                                 labelleft=False,
                                 right=False,
                                 labelright=False,
                                 top=False,
                                 labeltop=False)

        noTicks(chromPanel2)
        noTicks(chromPanel3)
        noTicks(chromPanel4)
        noTicks(chromPanel5)
        noTicks(chromPanel7)
        noTicks(chromPanel8)
        noTicks(chromPanel9)
        noTicks(chromPanel10)
        noTicks(chromPanel12)
        noTicks(chromPanel13)
        noTicks(chromPanel14)
        noTicks(chromPanel15)
        noTicks(chromPanel17)
        noTicks(chromPanel18)
        noTicks(chromPanel19)
        noTicks(chromPanel20)
        noTicks(chromPanel22)

        plotByChrom(vcf_dictionary, 'chr1', chromPanel1)
        plotByChrom(vcf_dictionary, 'chr2', chromPanel2)
        plotByChrom(vcf_dictionary, 'chr3', chromPanel3)
        plotByChrom(vcf_dictionary, 'chr4', chromPanel4)
        plotByChrom(vcf_dictionary, 'chr5', chromPanel5)
        plotByChrom(vcf_dictionary, 'chr6', chromPanel6)
        plotByChrom(vcf_dictionary, 'chr7', chromPanel7)
        plotByChrom(vcf_dictionary, 'chr8', chromPanel8)
        plotByChrom(vcf_dictionary, 'chr9', chromPanel9)
        plotByChrom(vcf_dictionary, 'chr10', chromPanel10)
        plotByChrom(vcf_dictionary, 'chr11', chromPanel11)
        plotByChrom(vcf_dictionary, 'chr12', chromPanel12)
        plotByChrom(vcf_dictionary, 'chr13', chromPanel13)
        plotByChrom(vcf_dictionary, 'chr14', chromPanel14)
        plotByChrom(vcf_dictionary, 'chr15', chromPanel15)
        plotByChrom(vcf_dictionary, 'chr16', chromPanel16)
        plotByChrom(vcf_dictionary, 'chr17', chromPanel17)
        plotByChrom(vcf_dictionary, 'chr18', chromPanel18)
        plotByChrom(vcf_dictionary, 'chr19', chromPanel19)
        plotByChrom(vcf_dictionary, 'chr20', chromPanel20)
        plotByChrom(vcf_dictionary, 'chr21', chromPanel21)
        plotByChrom(vcf_dictionary, 'chr22', chromPanel22)

    def histogramPanels(query_type, sample_name):
        panelWidth = 8
        panelHeight = 6

        # plot panels
        histPanel1 = plt.axes([0.05, 0.04, panelWidth / figureWidth, panelHeight / figureHeight])
        histPanel2 = plt.axes([0.53, 0.04, panelWidth / figureWidth, panelHeight / figureHeight])
        histPanel3 = plt.axes([0.05, 0.52, panelWidth / figureWidth, panelHeight / figureHeight])
        histPanel4 = plt.axes([0.53, 0.52, panelWidth / figureWidth, panelHeight / figureHeight])

        histPanel1.set_title('TP truth')
        histPanel2.set_title('TP query')
        histPanel3.set_title('FN')
        histPanel4.set_title('FP')

        histPanel1.set_xlabel('VAF')
        histPanel2.set_xlabel('VAF')
        histPanel3.set_xlabel('VAF')
        histPanel4.set_xlabel('VAF')
        histPanel1.set_ylabel('count')
        histPanel2.set_ylabel('count')
        histPanel3.set_ylabel('count')
        histPanel4.set_ylabel('count')

        plotHistograms(vcf_dictionary, histPanel1, histPanel2, histPanel3, histPanel4)

        if query_type == 'clairs':
            plt.suptitle('ClairS' + sample_name, fontsize=14)
        elif query_type == 'deepsomatic':
            plt.suptitle('DeepSomatic' + sample_name, fontsize=14)


    if plot_type == 'chromosome':
        chromPanels()
    elif plot_type == 'histogram':
        histogramPanels(query_type, sample_name)



    plt.savefig(out_file_name, dpi=600)
