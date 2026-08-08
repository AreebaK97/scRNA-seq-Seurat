## Dataset

Public GEO dataset:

- **GSM4138872**
- **GSM4138873**
- **GSM4138874**
- **GSM4138875**

### Dataset Note

The processed count matrices provided by the original study have already been filtered to remove mitochondrial and ribosomal genes.

Therefore, mitochondrial quality control cannot be calculated from these matrices.

## Preprocessing

The data were processed to prepare the Seurat objects for dimensionality reduction and downstream analysis.
The expression data was normalized to account for differences in sequencing depth between cells. Highly variable genes were also identified so that downstream analysis focuses on those genes that show the most variation between cells. These genes are more informative for distinguishing different cell types. 
The selected genes were centered and scaled so that genes are placed on a comparable scale before dimensionality reduction.

### Principal Component Analysis (PCA)
PCA was used to reduce the dimensionality of gene expression data. The number of principal components (PCs) used for downstream analysis was determined by examining the elbow plot. The plot showed that the variation started to decrease after PC 10. Therefore, the first 12 PCs were retained for further downstream analysis.

![PCA Elbow Plot](figures/preprocessing/BMMC_T1_elbow_plot.png)

The first 12 principal components were used for subsequent neighbourhood construction, clustering, UMAP visualization and doublet detection. 
