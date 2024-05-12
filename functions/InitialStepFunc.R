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
library(patchwork)
library(DT)
library(purrr)
library(pheatmap)
library(htmlwidgets)
library(networkD3)
library(edgeR)
library(stringr)
library(sctransform)
library(gridExtra)
library(grid)


normalizeData <- function(obj, method, using_log, norm_parameter = NULL, sct_parameter = NULL) {
    if (method == "sctransform") {
        if (using_log) {
            obj <- NormalizeData(obj, normalization_method = "LogNormalize")
        }
        obj <- SCTransform(obj, ncells = sct_parameter, vars.to.regress = "percent.mt", verbose = FALSE)
    } else {
        obj <- NormalizeData(obj, normalization_method = method, scale.factor = norm_parameter)
        obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    }
    all.genes <- rownames(obj)
    obj <- ScaleData(obj, features = all.genes)
    return(obj)
}

runPCA <- function(obj, num_pcs) {
    num_pcs <- min(num_pcs, ncol(obj))
    if ("sctransform" %in% names(obj@assays)) {
        obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = num_pcs)
    } else {
        obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = num_pcs, verbose = FALSE)
    }
    return(obj)
}

findVariableFeatures <- function(obj, num_features) {
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = num_features)
    top10 <- head(VariableFeatures(obj), 10)
    return(list(object = obj, top10 = top10))
}

plotVariableFeatures <- function(obj, top10) {
    plot <- LabelPoints(plot = VariableFeaturePlot(obj), points = top10, repel = TRUE)
    return(plot)
}
