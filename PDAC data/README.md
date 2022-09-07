Single-cell analysis of the pancreatic cancer microbiome 
================
Bassel Ghaddar
7/31/2022

### Introduction

We used SAHMI to analyze tumor-microbiome interactions in two pancreatic cancer cohorts: scPDA1 (Peng et al., 2019) and scPDA2 (Steele et al., 2020). We also compare metagnomic read capture in single-cell RNA-seq, whole genome-seq, 16S-rDNA-seq, and bulk RNA-seq.  

File descriptions: 
- pdac_sahmi.RDS: SAHMI analysis results for scPDA1 and scPDA2, including metagenomic read counts and k-mer correlation test and cell-line quantile test results. 

- rutgers_sahmi.RDS: metagenomic read counts for all in-house samples analyzed by 16S-rDNA-seq, single-cell RNA-seq, bulk RNA-seq, and whole genome-seq. Raw sequencing files can be found at phs003035.v1.p1. 

- peng_microbiome.RDS: processed single-barcode by taxa counts matrix for scPDA1. 

- steele_microbiome.RDS: processed single-barcode by taxa counts matrix for scPDA2. 

Please contact Bassel Ghaddar (bassel.ghaddar@gmail.com) or Subhajyoti De (subhajyoti.de@rutgers.edu) for any questions. 

### References
Coming soon.

Ghaddar, B., Blaser, M.J., and De, S. (2022). Denoising sparse microbial signals from single-cell sequencing of mammalian host tissues. BioRxiv. doi: 10.1101/2022.06.29.498176

Peng, J., Sun, B.F., Chen, C.Y., Zhou, J.Y., Chen, Y.S., Chen, H., Liu, L., Huang, D., Jiang, J., Cui, G.S., et al. (2019). Single-cell RNA-seq highlights intra-tumoral heterogeneity and malignant progression in pancreatic ductal adenocarcinoma. Nat. Cell Res.

Steele, N.G., Carpenter, E.S., Kemp, S.B., Sirihorachai, V.R., The, S., Delrosario, L., Lazarus, J., Amir, E. ad D., Gunchick, V., Espinoza, C., et al. (2020). Multimodal mapping of the tumor and peripheral blood immune landscape in human pancreatic cancer. Nat. Cancer 1, 1097–1112.

