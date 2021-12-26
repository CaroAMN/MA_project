#!/bin/bash
# ----------------
# |   Purpose    |
# ----------------
# Cluster maf files into a gbm, primary and recurrent 

for file in vcf_filtered_maf/*.maf; do
  filename="$(basename -- $file)"
  awk 'NR>2' $file > vcf_filtered_maftools/"$filename"
done
# Cat vcfs to use them into maftools
# dont forget to add a header!
cat vcf_filtered_maftools/* > all_vcfs_filtered.maf
find vcf_filtered_maftools/* -name '*_T_*' -exec cat {} + > primary_vcfs_filtered.maf
find vcf_filtered_maftools/* -name '*_R_*' -exec cat {} + > recurrent_vcfs_filtered.maf
