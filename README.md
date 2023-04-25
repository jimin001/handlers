# handlers

## Scripts for handling BAM, FASTA and VCF files.

#### plot_vcf.py
Plots VAF from VCFs. 2 plot style options: scatterplot by chromosome position colored by FN, FP, TP, or histogram with 4 panels for FN, FP, TP(truth), TP(query)

run command:
```sh
  python3 plot_vcf.py \
  --queryvcf query.vcf.gz \
  --queryvcfindex query.vcf.gz.tbi \
  --truthvcf truth.vcf.gz \
  --truthvcfindex truth.vcf.gz.tbi \
  --features sompy.features.csv \
  --querytype 'variant caller' \
  --plot 'plot type' \
  --samplename 'name of sample' \
  -o output_file_name.png
  ```
