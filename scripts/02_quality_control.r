###########################################################
# Script : 01_quality_control.r
# Description : 
# Load seurat objects and perform QC
###########################################################

library(here)
library(Seurat)
library(readr)

source(here("scripts", "config.r"))

sample_metadata <- read_tsv(here("metadata", "sample_metadata.tsv"))   

quality_control <- function(srt_obj, sample){
    srt_obj[["percent.mt"]] <- PercentageFeatureSet(
        srt_obj, pattern = "^MT-"
    )

    # head(srt_obj@meta.data)

    vln_plt <- VlnPlot(srt_obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

    ggsave(filename = here("figures", "qc",
        paste0(sample, "_qc.png")), plot = vln_plt)

    scatter_plt1 <- FeatureScatter(srt_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    scatter_plt2 <- FeatureScatter(srt_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")

    ggsave(filename = here("figures", "qc",
        paste0(sample, "_scatter_plt_CountvsFeature.png")), plot = scatter_plt1)

    ggsave(filename = here("figures", "qc",
        paste0(sample, "_scatter_plt_CountvsMT.png")), plot = scatter_plt2)


    srt_obj <- subset(srt_obj, subset = nFeature_RNA > min_feature & nFeature_RNA < max_feature & percent.mt < 5)

    return(srt_obj)
}

for (i in 1:nrow(sample_metadata)){
    sample_obj <- readRDS(here("data", "processed", paste0(sample_metadata$sample[i], ".rds")))
    # print(sample_obj)
    sample_obj <- quality_control(sample_obj, sample_metadata$sample[i])
    saveRDS(sample_obj,
        here("data", "processed", paste0(sample_metadata$sample[i], "_QC.rds")))
}

sample_obj <- readRDS(here("data", "processed", "BMMC_T1.rds"))

sample_obj[["percent.mt"]] <- PercentageFeatureSet(
        sample_obj, pattern = "^MT-"
    )
head(sample_obj@meta.data)
summary(sample_obj$percent.mt)
