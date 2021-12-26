#!/bin/bash
# ----------------
# |   Purpose    |
# ----------------
# Convert vcf to maf files using vcf2maf 

for file in vcf_filtered/*.vcf; do
	filename="$(basename -- $file)"
	filename="${filename%.*}"
	tumorID=$(basename $file | cut -d'_' -f 2)
  normalID=$(basename $file | cut -d'_' -f 5)
	perl vcf2maf.pl --input-vcf "$file" \
	--output-maf vcf_filtered_maf/"$filename".maf \
	--tumor-id "$tumorID" \
	--normal-id "$normalID" \
	--vcf-tumor-id TUMOR \
	--vcf-normal-id NORMAL \
	--vep-path ~/opt/anaconda3/envs/vcf2maf/bin
done
