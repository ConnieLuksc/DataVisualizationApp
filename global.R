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
library(patchwork)
library(DT)
library(purrr)
library(pheatmap)

library(networkD3)
library(edgeR)


# Read in file and perform validation.
load_seurat_obj <- function(path) {
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
  if (!inherits(obj, "Seurat")) {
    errors <- c(errors, "File is not a seurat object")
    return(errors)
  }

  # mock annotation column in metadata
  annotation <- list(ambiguous = list(0, 8), cd45.1 = list(1, 3), cd45.2 = list(2, 4, 6), cd45.3 = list(5, 7))
  new_metadata <- rep(NA, nrow(obj))
  # assign cluster with annotation
  for (key in names(annotation)) {
    clusters <- unlist(annotation[[key]])
    for (cluster in clusters) {
      new_metadata[obj$seurat_clusters == cluster] <- key
    }
  }
  obj <- AddMetaData(obj, metadata = new_metadata, col.name = "Annotation")

  return(obj)
}


create_metadata_UMAP <- function(obj, col, pc, resolution, values) {
  tryCatch(
  {
    obj <- FindNeighbors(obj, dims = 1:pc)
    obj <- FindClusters(obj, resolution = resolution)
    obj <- RunUMAP(obj, dims = 1:pc)
    umap <- DimPlot(obj, pt.size = .1, label = FALSE, label.size = 4, group.by = col, reduction = "umap") # nolint: line_length_linter.
    remove_modal_spinner()
    values$obj <- obj
    values$umap <- umap
    return(umap)
  },
    error = function(err) {
      umap <- ggplot() +
        theme_void() +
        geom_text(aes(x = 0.5, y = 0.5, label = str_wrap(err$message, width = 20)), size = 12, color = "gray73", fontface = "bold") + # nolint: line_length_linter.
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
      values$obj <- NULL
      remove_modal_spinner()
      return(umap)
    }
  )
}

create_feature_plot <- function(obj, genes, values) {
    FP <- NULL
    if(!is.null(genes)){
        genes_list <- strsplit(genes, "\\s+")
        obj <- AddModuleScore(obj, features = genes_list, name = "module")
        FP <- FeaturePlot(obj, features = "module1", label = TRUE, repel = TRUE)   
    }


    # if (gene %in% rownames(obj)) {
        # FP <- FeaturePlot(obj, features = gene, pt.size = 0.001, combine = FALSE)
    # } else {
    #     FP <- ggplot() +
    #     theme_void() +
    #     geom_text(aes(x = 0.5, y = 0.5, label = str_wrap("Gene doesn't exist", width=20)), size = 12, color = "gray73", fontface = "bold") +
    #     theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    # }
    # values$feature <- FP
    return(FP)
}

create_violin_plot <- function(obj, genes, values, ncol, pt.size) {
    combined_plot <- NULL
    VP <- list()

    if(!is.null(genes)){
        genes_list <- strsplit(genes, "\\s+")
        for (gene in genes_list) {
            VP[[gene]] <- VlnPlot(obj, features = gene, combine = TRUE)
        }
        combined_plot <- cowplot::plot_grid(plotlist = VP, ncol = length(genes_list))
    }    

    # if (gene %in% rownames(obj)) {
    #     VP <- VlnPlot(obj, features = gene, pt.size = 0, combine = FALSE)
    # }else{
    #     VP <- ggplot() +
    #     theme_void() +
    #     geom_text(aes(x = 0.5, y = 0.5, label = str_wrap("Gene doesn't exist", width=20)), size = 12, color = "gray73", fontface = "bold") +
    #     theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    # }
    # values$violin <- VP
    return(combined_plot)
}

create_mds_plot <- function(obj, values) {
  # aggregate
  all_genes <- AggregateExpression(obj)

  # create a deglist object
  gene_names <- rownames(obj@assays$RNA)
  dge_data <- DGEList(counts = all_genes$RNA)
  dge_data$genes <- gene_names

  # create MDS plot
  set.seed(123)
  colors <- rainbow(length(colnames(all_genes$RNA)))
  # mdsPlot <- plotMDS(dge_data, col = colors)
  mdsPlot <- plotMDS(dge_data, col = colors, pch=20, cex=2)
  legend("topright", legend = colnames(all_genes$RNA), col = colors, pch = 20, cex = 0.8, pt.cex = 0.8, title = "Group")
  title("MDS Plot")
  values$mds <- mdsPlot
  return (mdsPlot)
}

# visualize annotation
create_annotation_UMAP <- function(obj, col, pc, resolution, values, annotation) {
  tryCatch(
  {
    obj <- FindNeighbors(obj, dims = 1:pc)
    obj <- FindClusters(obj, resolution = resolution)
    obj <- RunUMAP(obj, dims = 1:pc)
    umap <- DimPlot(obj,
                    pt.size = .1, label = FALSE,
                    label.size = 4, group.by = "Annotation",
                    reduction = "umap"
    )
    remove_modal_spinner()
    values$obj <- obj
    values$umap_annotation <- umap
    return(umap)
  },
    error = function(err) {
      umap <- ggplot() +
        theme_void() +
        geom_text(
          aes(
            x = 0.5, y = 0.5,
            label = str_wrap(err$message, width = 20)
          ),
          size = 12,
          color = "gray73", fontface = "bold"
        ) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
      values$obj <- NULL
      # remove_modal_spinner()
      return(umap)
    }
  )
}

# sankey plot
create_sankey_plot <- function(obj, values, pc, resolution, annotation, cluster_num) {
  # obj <- FindNeighbors(obj, dims = 1:pc)
  # obj <- FindClusters(obj, resolution = resolution)
  cluster_info <- as.data.frame(obj$seurat_clusters)
  annotation_info <- as.data.frame(obj$Annotation)
  transitions <- data.frame(From = character(), To = character())

  # record transition
  for (i in 0:nrow(cluster_info)) {
    cluster <- cluster_info[i, 1]
    annotation <- annotation_info[i, 1]
    transition <- data.frame(Cluster = cluster, Annotation = annotation)
    transitions <- rbind(transitions, transition)
  }

  # calculate frequency of transition
  transition_counts <- table(transitions)

  # create transition
  transition_data <- data.frame(transition_counts)
  colnames(transition_data) <- c("From", "To", "Value")
  transition_data$From <- as.integer(transition_data$From) - 1
  transition_data$To <- as.character(transition_data$To)
  for (i in seq_along(transition_data$To)) {
    transition_data$To[i] <- as.character(which(colnames(transition_counts) == transition_data$To[i]))
  }
  transition_data$To <- as.integer(transition_data$To) + as.integer(length(rownames(transition_counts))) - 1

  # Nodes parameter
  all_nodes <- data.frame(name = c(as.character(rownames(transition_counts)), as.character(colnames(transition_counts))))

  # create sankey
  sankey <- sankeyNetwork(
    Links = transition_data,
    Source = "From",
    Nodes = all_nodes,
    Target = "To",
    Value = "Value",
    units = "Cell Counts",
    width = 600,
    height = 400
  )
  values$sankey <- sankey
  return(sankey)

}
