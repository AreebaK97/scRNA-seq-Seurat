##################################################################
# Script : 04_doublet_removal.r
# Description : Identify and remove doublets using DoubletFinder
# Neighbour graph and cluster construction
# Optimal pK value determination
##################################################################


library(here)
library(Seurat)
library(readr)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
source(here("scripts", "config.r"))

sample_metadata <- read_tsv(here("metadata", "sample_metadata.tsv"))

doublet_removal <- function(srt_obj, sample_name){
    srt_obj <- FindNeighbors(
        srt_obj,
        dims = pcs_to_use
    )

    srt_obj <- FindClusters(
        srt_obj,
        resolution = cluster_resolution
    )

    srt_obj <- RunUMAP(
        srt_obj,
        dims = pcs_to_use
    )

    umap_before <- DimPlot(
        srt_obj,
        reduction = "umap",
        label = TRUE
    )
    ggsave(
        filename = here(
            "figures",
            "doublets",
            paste0(sample_name, "_umap_before_doublet_removal.png")
        ),
        plot = umap_before,
        width = 8,
        height = 6,
        dpi = 300,
        bg = "white")

    
    pk_value <- find_pk_value(srt_obj, sample_name)

    annotations <- srt_obj$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)

    nExp_poi <- round(doublet_rate * ncol(srt_obj))
    nExp_poi.adj <- round(
        nExp_poi * (1 - homotypic.prop)
    )


    srt_obj <- doubletFinder(
        srt_obj,
        PCs = pcs_to_use,
        pN = 0.25,
        pK = pk_value,
        nExp = nExp_poi.adj,
        reuse.pANN = NULL,
        sct = FALSE
    )

    saveRDS(srt_obj, here("data", "processed", paste0(sample_name, "_doubletfinder_results.rds")))

    classification_col <- grep(
        "DF.classifications",
        colnames(srt_obj@meta.data),
        value = TRUE
    )

    doublet_plt <- DimPlot(
        srt_obj,
        reduction = "umap",
        group.by = classification_col
    )

    ggsave(
        filename = here(
            "figures",
            "doublets",
            paste0(sample_name, "_doublet_classification.png")
        ),
        plot = doublet_plt,
        width = 8,
        height = 6,
        dpi = 300,
        bg = "white")


    singlet_cells <- rownames(
        srt_obj@meta.data[
            srt_obj@meta.data[[classification_col]] == "Singlet",
        ]
    )

    srt_obj_clean <- subset(
        srt_obj,
        cells = singlet_cells
    )


    umap_after <- DimPlot(
        srt_obj_clean,
        reduction = "umap",
        label = TRUE
    )

    ggsave(
        filename = here(
            "figures",
            "doublets",
            paste0(sample_name, "_umap_after_doublet_removal.png")
        ),
        plot = umap_after,
        bg = "white",
        width = 8,
        height = 6,
        dpi = 300
)
    return(srt_obj_clean)
}

find_pk_value <- function(srt_obj, sample){
    sweep_results <- paramSweep(
        srt_obj,
        PCs = pcs_to_use,
        sct = FALSE
    )

    sweep_stats <- summarizeSweep(
        sweep_results,
        GT = FALSE
    )

    bcmvn <- find.pK(sweep_stats)
    pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
    pK <- as.numeric(as.character(pK))

    pk_plot <- ggplot(
        bcmvn,
        aes(
            x = as.numeric(as.character(pK)),
            y = BCmetric
        )
    ) +
        geom_point() +
        geom_line() +
        theme_classic()

    ggsave(
        here(
            "figures",
            "doublets",
            paste0(sample, "_pK_selection.png")
        ),
        pk_plot,
        bg = "white"
    )

    return(pK)

}



for (i in 1:nrow(sample_metadata)){ 
    sample_obj <- readRDS(here("data", "processed", paste0(sample_metadata$sample[i], "_preprocessed.rds")))
    # print(sample_obj)
    sample_obj <- doublet_removal(sample_obj, sample_metadata$sample[i])
    saveRDS(sample_obj,
        here("data", "processed", paste0(sample_metadata$sample[i], "_doublets_removed.rds")))
}

