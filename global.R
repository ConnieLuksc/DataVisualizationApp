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
library(stringr)
library(sctransform)

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

  return(obj)
}


create_metadata_UMAP <- function(obj, pc, resolution, values) {
  tryCatch(
  {
    obj <- FindNeighbors(obj, dims = 1:pc)
    obj <- FindClusters(obj, resolution = resolution)
    obj <- RunUMAP(obj, dims = 1:pc)
    colors <- rainbow(length(levels(Idents(obj))))
    umap <- DimPlot(obj, pt.size = .1, label = FALSE, label.size = 4, reduction = "umap", cols=colors) 
    # umap <- DimPlot(obj, pt.size = .1, label = FALSE, label.size = 4, group.by = col, reduction = "umap", cols=colors) 
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
        combined_genes <- paste(genes, collapse = " ")
        wrapped_title <- str_wrap(combined_genes, width = 10)
        obj <- AddModuleScore(obj, features = genes_list, name = combined_genes)
        FP <- FeaturePlot(obj, features = paste0(combined_genes, "1"), label = TRUE, repel = TRUE)+labs(title = wrapped_title) +
      theme(plot.title = element_text(size = 12)) 
    }

    return(FP)
}

create_violin_plot <- function(obj, genes, values, ncol, pt.size) {
    VP <- NULL

    if(!is.null(genes)){
        genes_list <- strsplit(genes, "\\s+")
        combined_genes <- paste(genes, collapse = " ")
        wrapped_title <- str_wrap(combined_genes, width = 10)
        obj <- AddModuleScore(obj, features = genes_list, name = combined_genes)
        VP <- VlnPlot(obj, features = paste0(combined_genes, "1"), pt.size = 0, combine = TRUE)+labs(title = wrapped_title) +
      theme(plot.title = element_text(size = 12))   
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
  legend("topright", legend = colnames(all_genes$RNA), col = colors, pch = 20, cex = 0.8, pt.cex = 0.8, title = "Group")
  title("MDS Plot")
  values$mds <- mdsPlot
  return (mdsPlot)
}

# visualize annotation
create_annotation_UMAP <- function(obj, pc, resolution, values, annotation) {
  tryCatch(
  {
    obj1 <- RunUMAP(obj, dims = 1:pc)
    umap <- DimPlot(obj1,
                    pt.size = .1, label = FALSE,
                    label.size = 4, group.by = annotation,
                    reduction = "umap"
    )
    remove_modal_spinner()
    # values$obj <- obj
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
      # remove_modal_spinner()
      return(umap)
    }
  )
}

# sankey plot
create_sankey_plot <- function(obj, values, pc, resolution, annotation, cluster_num) {
  cluster_info <- as.data.frame(obj$seurat_clusters)
  # annotation_info <- as.data.frame(obj$Annotation)
  annotation_info <- as.data.frame(obj[[annotation]])
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
