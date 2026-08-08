###########################################################
# Script : 03_preprocessing.r
# Description : 
# Normalize scRNA-seq data, identify highly variable genes,
# scale the data, and perform PCA.
###########################################################


library(here)
library(Seurat)
library(readr)
library(ggplot2)

source(here("scripts", "config.r"))

sample_metadata <- read_tsv(here("metadata", "sample_metadata.tsv"))


preprocessing <- function(srt_obj, sample_name){
    
    ## Normalization
    srt_obj <- NormalizeData(srt_obj, normalization.method = "LogNormalize", scale.factor = 1000)

    ## Feature Selection and Visualize
    srt_obj <- FindVariableFeatures(srt_obj, selection.method = "vst", nfeatures = n_variable_features)
    top10_features <- head(VariableFeatures(srt_obj), 10)

    feature_plot <- VariableFeaturePlot(srt_obj)
    feature_plot <- LabelPoints(plot = feature_plot, points = top10_features, repel = TRUE)

    ggsave(
        filename = here(
            "figures",
            "preprocessing",
            paste0(sample_name, "_variable_features.png")
        ),
        plot = feature_plot,
        width = 8,
        height = 6,
        dpi = 300,
        bg = "white")

    ## Scale Data
    srt_obj <- ScaleData(srt_obj, features = rownames(srt_obj))

    ## Run PCA and Visualize
    srt_obj <- RunPCA(srt_obj, features = VariableFeatures(srt_obj), npcs = n_pcs)

    # Elbow Plot 
    pca_elbow_plt <- ElbowPlot(srt_obj)
    ggsave(
        filename = here(
            "figures",
            "preprocessing",
            paste0(sample_name, "_pca_elbow_plot.png")
        ),
        plot = pca_elbow_plt,
        width = 8,
        height = 6,
        dpi = 300,
        bg = "white")

    # PCA Loadings Plot
    pca_loadings_plt <- VizDimLoadings(srt_obj, dims = 1:2, reduction = "pca")
    ggsave(
        filename = here(
            "figures",
            "preprocessing",
            paste0(sample_name, "_pca_loadings.png")
        ),
        plot = pca_loadings_plt, 
        width = 8,
        height = 6,
        dpi = 300,
        bg = "white")
    

    # PCA
    dim_plt <- DimPlot(srt_obj, reduction = "pca") + NoLegend()
    ggsave(
        filename = here(
            "figures",
            "preprocessing",
            paste0(sample_name, "_pca_dim_plot.png")
        ),
        plot = dim_plt,
        width = 8,
        height = 6,
        dpi = 300,
        bg = "white")
    
    
   png(
    filename = here(
        "figures",
        "preprocessing",
        paste0(sample_name, "_pca_heatmap.png")
    ),
        width = 2100,
        height = 2500,
        res = 300,
        bg = "white"
    )

    DimHeatmap(
        srt_obj,
        dims = 1:12,
        cells = 500,
        balanced = TRUE
    )

    dev.off()

    return(srt_obj)
}


for (i in 1:nrow(sample_metadata)){
    sample_obj <- readRDS(here("data", "processed", paste0(sample_metadata$sample[i], "_QC.rds")))
    # print(sample_obj)
    sample_obj <- preprocessing(sample_obj, sample_metadata$sample[i])
    saveRDS(sample_obj,
        here("data", "processed", paste0(sample_metadata$sample[i], "_preprocessed.rds")))
}


