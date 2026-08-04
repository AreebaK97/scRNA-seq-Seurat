###########################################################
# Script : 01_load_data.r
# Description : 
# Load raw RNA-seq Data and create seurat objects
###########################################################

library(here)
library(Seurat)
library(tidyverse)
library(dplyr)
library(readr)

source(here("scripts", "config.r"))

sample_metadata <- read_tsv(
    here("metadata", "sample_metadata.tsv")
)

# print(sample_metadata)

##### Seurat Objects #####

create_seurat_obj <- function(file, sample, donor, replicate, sex){

    counts <- readRDS(here("data", "raw", file))

    obj <- CreateSeuratObject(counts, project = sample,
        min.cells = 3, min.features = min_feature)

    obj$donor <- donor
    obj$sex <- sex
    obj$replicate <- replicate

    return(obj)
}

for (i in 1:nrow(sample_metadata)){
    seurat_obj <- create_seurat_obj(sample_metadata$file[i],
        sample_metadata$sample[i],
        sample_metadata$donor[i],
        sample_metadata$replicate[i],
        sample_metadata$sex[i]
    )
    # print(seurat_obj)

    saveRDS(seurat_obj,
        here("data", "processed", paste0(sample_metadata$sample[i], ".rds")))
}
