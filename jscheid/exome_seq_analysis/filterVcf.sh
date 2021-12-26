#!/bin/bash
# ----------------
# |   Purpose    |
# ----------------
# Use bcftools to filter the each vcf file

for file in vcf/*.vcf; do
  filename="$(basename -- $file)"
  bcftools view -i 'FILTER="PASS"' $file > vcf_filtered/filtered"$filename"
done
