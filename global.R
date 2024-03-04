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
library(htmlwidgets)

library(networkD3)
library(edgeR)

# set.seed(123)
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
    colors <- rainbow(length(levels(Idents(obj))))
    umap <- DimPlot(obj, pt.size = .1, label = FALSE, label.size = 4, group.by = col, reduction = "umap", cols=colors) # nolint: line_length_linter.
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
        obj <- AddModuleScore(obj, features = genes_list, name = "feature_plot")
        # colors <- rainbow(length(levels(Idents(obj))))
        FP <- FeaturePlot(obj, features = "feature_plot1", label = TRUE, repel = TRUE)   
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

create_violin_plot <- function(obj, genes = NULL, values, ncol = NULL, pt.size) {
    VP <- NULL
    colors <- rainbow(length(levels(Idents(obj))))

    if(!is.null(genes) && length(genes) > 0){
        genes_list <- strsplit(genes, "\\s+")
        obj <- AddModuleScore(obj, features = genes_list, name = "violin_plot")
        VP <- VlnPlot(obj, features = "violin_plot1", pt.size = pt.size, combine = TRUE, cols = colors)
    } else {
        # Case when genes input is not provided
        VP <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = ncol, pt.size = pt.size, cols = colors)
    }

    return(VP)
}

create_mds_plot <- function(obj, values) {
  # aggregate
  all_genes <- AggregateExpression(obj)

  # create a deglist object
  gene_names <- rownames(obj@assays$RNA)
  dge_data <- DGEList(counts = all_genes$RNA)
  dge_data$genes <- gene_names

  # create MDS plot
  # set.seed(123)
  colors <- rainbow(length(colnames(all_genes$RNA)))
  values$color <- colors
  # mdsPlot <- plotMDS(dge_data, col = colors)
  mdsPlot <- plotMDS(dge_data, col = colors, pch=20, cex=2)
  # legend("topright", legend = colnames(all_genes$RNA), col = colors, pch = 20, cex = 0.8, pt.cex = 0.8, title = "Group")
  legend("topright", legend = gsub("^g", "", colnames(all_genes$RNA)), col = colors, pch = 20, cex = 0.8, pt.cex = 0.8, title = "Group")
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

  unique_from <- unique(transition_data$From)
  # colors <- rainbow(length(unique_from))
  colors <- values$color
  names(colors) <- unique_from
  print(colors)

  transition_data$Color <- colors[as.character(transition_data$From)]

  # Nodes parameter
  all_nodes <- data.frame(name = c(as.character(rownames(transition_counts)), as.character(colnames(transition_counts))))

  # Change the color of nodes
  all_nodes$group <- as.character(all_nodes$name)
  for (node in unique_from) {
    if (node %in% all_nodes$name) {
      all_nodes$color[all_nodes$name == node] <- colors[node+1]
    }
  }
  colors_js_array <- paste0("['", paste(values$color, collapse = "','"), "']")
  colourScale_js <- sprintf("d3.scaleOrdinal().domain(d3.range(1,%s + 1)).range(%s)", length(values$color), colors_js_array)

  # create sankey
  sankey <- sankeyNetwork(
    Links = transition_data,
    Source = "From",
    Nodes = all_nodes,
    Target = "To",
    Value = "Value",
    units = "Cell Counts",
    NodeID = "name",
    NodeGroup = "color",
    colourScale = JS(colourScale_js),
    width = 600,
    height = 400
  )

  sankey <- onRender(
  sankey,
  '
  function(el, x) {
    d3.select(el).selectAll(".link")
      .style("stroke", function(d) { return d.source.color; });
  }
  '
)

  values$sankey <- sankey
  return(sankey)

}
