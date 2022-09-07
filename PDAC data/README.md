Single-cell analysis of the pancreatic cancer microbiome 
================
Bassel Ghaddar
7/31/2022

### Introduction

We used SAHMI to analyze tumor-microbiome interactions in pancreatic cancer cohorts and report the processed metagenomic counts data for all datasets analyzed in Ghaddar et al., 2022 in review. 

File descriptions: 
- rutgers_sahmi.RDS: SAHMI metagenomic read count reports for all in-house samples analyzed by 16S-rDNA-seq, single-cell RNA-seq, bulk RNA-seq, and whole genome-seq. Raw sequencing files can be found at phs003035.v1.p1. 

- pdac_sahmi.RDS: SAHMI metagenomic read count reports for single-cell sequencing data from two cohorts (scPDA1 and scPDA2). 

- peng_microbiome.RDS: processed single-barcode by taxa counts matrix for scPDA1. 

- steele_microbiome.RDS: processed single-barcode by taxa counts matrix for scPDA2. 

Please contact Bassel Ghaddar (bassel.ghaddar@gmail.com) or Subhajyoti De (subhajyoti.de@rutgers.edu) for any questions. 

### References

Ghaddar, B., Blaser, M.J., and De, S. (2022). Denoising sparse microbial signals from single-cell sequencing of mammalian host tissues. BioRxiv. doi: 10.1101/2022.06.29.498176

Peng, J., Sun, B.F., Chen, C.Y., Zhou, J.Y., Chen, Y.S., Chen, H., Liu, L., Huang, D., Jiang, J., Cui, G.S., et al. (2019). Single-cell RNA-seq highlights intra-tumoral heterogeneity and malignant progression in pancreatic ductal adenocarcinoma. Nat. Cell Res.

Steele, N.G., Carpenter, E.S., Kemp, S.B., Sirihorachai, V.R., The, S., Delrosario, L., Lazarus, J., Amir, E. ad D., Gunchick, V., Espinoza, C., et al. (2020). Multimodal mapping of the tumor and peripheral blood immune landscape in human pancreatic cancer. Nat. Cancer 1, 1097–1112.

