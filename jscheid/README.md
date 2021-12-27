# Thesis scripts
**Title:** Integration of multi-omics data to investigate differences and overlaps of primary and recurrent glioblastoma

## Table of contents
* [Deconvolution analysis](#deconvolution-analysis)
* [DNA methylation analysis](#dna-methylation-analysis)
* [Exome sequencing analysis](#exome-sequencing-analysis)
* [HLA ligandomics analysis](#hla-ligandomics-analysis)
* [RNA sequencing analysis](#rna-sequencing-analysis)

## Deconvolution analysis
Deconvolution of bulk RNA-seq data using published celltype transcript profiles of [Neftel et al. 2019](https://www.sciencedirect.com/science/article/pii/S0092867419306877#) in [CIBERSORTx](https://cibersortx.stanford.edu) and [Scaden](https://scaden.ims.bio).
The output of both tools was processed using `neftelPostProcessing.R` and visualized with `deconvolutionVisualizations.ipynb` to infer the cell composition in the bulk RNA-seq data.

## DNA methylation analysis
`SummaryCNVs.R` reads in `.igv` files containing information about the log2-transformed copy-number ratio. The copy-number ratio is measured by the sum of methylated and unmethylated signal intensities of the malignant samples against an [in-house database](https://pubmed.ncbi.nlm.nih.gov/29967940/) of healthy reference cohort. The `.igv` files were summarized to visualize the DNA methylation profile of primary and recurrent GBM cohorts. Additionaly, a horizon plot visualizes the per-sample DNA methylation profile.

## Exome sequencing analysis
`.vcf` files were filtered using `filterVcf.sh` and converted into `.maf` files using `vcf2maf.sh`. The resulting `.maf` files were combined into 3 `.maf` databases: Primary, recurrent and the whole GBM cohort using `clusterVcf.sh`. Those 3 databases served as input for `maftools.R`, visualizes produces oncoplots, somatic interactions, oncogeny pathways etc.

## HLA ligandomics analysis
`ligandomicsAnalysis.R` reads in an excel file containing HLA class 1 or class 2 peptides with their respective sample identifier. For HLA class 1 and class 2 it plots:
* Peptide yield
* Length distribution 
* Waterfall plots against a benign database
* Saturation

## RNA sequencing analysis
`rawCountToTPM.R` was used to transform the RNA-seq data into transcripts-per-million (TPM) in order to plug the data into `CIBERSORTx` and `Scaden`.
`RNASeqAnalysis.R` reads in rlog or vst transformed count data as well as differentially expressed genes returned by [DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). D.e. genes were clustered in a heatmap. Additionally, a PCA was performed to investigate potential batches in the data.
