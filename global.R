# library(gdata)
library(tools)
library(Seurat)
library(dplyr)
library(shinyjs)
library(shiny)
library(DT)
library(shinydashboard)
library(shinydashboardPlus)
library(ggplot2)
library(shinybusy)
library(glue)
library(markdown)
library(ggthemes)
library(stringr)


# Read in file and perform validation.
load_seurat_obj <- function(path){
    errors <- c()
    # check file extension
    if (!tolower(tools::file_ext(path)) == "rds") { # ignores case
        errors <- c(errors, "Invalid rds file.")
        return(errors)
    }

    # try to read in file
    tryCatch(
        {
        obj <- readRDS(path)
        },
        error = function(e) {
            errors <- c(errors, "Invalid rds file.")
            return(errors)
        }
    )

    # Validate obj is a seurat object
    if (!inherits(obj, "Seurat")){
        errors <- c(errors, "File is not a seurat object")
        return(errors)
    }

    return(obj)
}


create_metadata_UMAP <- function(obj, col, pc, resolution, values){
    tryCatch({
        obj <- FindNeighbors(obj, dims = 1:pc)
        obj <- FindClusters(obj, resolution = resolution)
        obj <- RunUMAP(obj, dims = 1:pc)
        umap <- DimPlot(obj, pt.size = .1, label = FALSE, label.size = 4, group.by = col, reduction = "umap")
        remove_modal_spinner()
        values$obj <- obj
        values$umap <- umap
        return(umap)
    }, error = function(err) {
        umap <- ggplot() +
            theme_void() +
            geom_text(aes(x = 0.5, y = 0.5, label = str_wrap(err$message, width=20)), size = 12, color = "gray73", fontface = "bold") +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
            values$obj <- NULL
            remove_modal_spinner()
        return(umap)
    })
}

create_feature_plot <- function(obj, gene, values) {
    if (gene %in% rownames(obj)) {
        FP <- FeaturePlot(obj, features = gene, pt.size = 0.001, combine = FALSE)
    } else {
        FP <- ggplot() +
        theme_void() +
        geom_text(aes(x = 0.5, y = 0.5, label = "Gene doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    }
    values$feature <- FP
    return(FP)
}

create_violin_plot <- function(obj, gene, values) {
    if (gene %in% rownames(obj)) {
        VP <- VlnPlot(obj, features = gene, pt.size = 0.001, combine = FALSE)
    }
    values$violin <- VP
    return(VP)
}