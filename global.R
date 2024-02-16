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
library(networkD3)


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


create_metadata_UMAP <- function(obj, col, pc, resolution, values) {
  tryCatch(
  {
    obj <- FindNeighbors(obj, dims = 1:pc)
    obj <- FindClusters(obj, resolution = resolution)
    obj <- RunUMAP(obj, dims = 1:pc)
    umap <- DimPlot(obj, pt.size = .1, label = FALSE, label.size = 4, group.by = col, reduction = "umap") # nolint: line_length_linter.
    # remove_modal_spinner()
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
      # remove_modal_spinner()
      return(umap)
    }
  )
}

create_feature_plot <- function(obj, gene, values) {
  if (gene %in% rownames(obj)) {
    FP <- FeaturePlot(obj, features = gene, pt.size = 0.001, combine = FALSE)
  } else {
    FP <- ggplot() +
      theme_void() +
      geom_text(aes(x = 0.5, y = 0.5, label = str_wrap("Gene doesn't exist", width = 20)), size = 12, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  values$feature <- FP
  return(FP)
}

create_violin_plot <- function(obj, gene, values, ncol, pt.size) {
  VP <- NULL
  if (gene %in% rownames(obj)) {
    VP <- VlnPlot(obj, features = gene, pt.size = 0.001, combine = FALSE)
  } else {
    VP <- ggplot() +
      theme_void() +
      geom_text(aes(x = 0.5, y = 0.5, label = str_wrap("Gene doesn't exist", width = 20)), size = 12, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  values$violin <- VP
  return(VP)
}

create_mds_plot <- function(obj, values) {
  # representative_cells <- data.frame()
  # for (cluster_id in unique(Idents(values$obj))) {
  #   cluster_cells <- Idents(obj) == cluster_id
  #   # print(cluster_cells)
  #   cluster_data <- obj[, cluster_cells]
  #   print(str(cluster_data))
  #   representative_cell <- as.data.frame(t(as.matrix(GetAssayData(object = cluster_data, assay = "RNA"))[1 ,])) # 选择第一个细胞作为代表性细胞，你也可以根据自己的需求选择
  #   print(as.character(rownames(representative_cell)))
  #   rownames(representative_cell) <- rownames(cluster_data)[1]
  #   # print(str(obj))
  #   representative_cells <- rbind(representative_cells, representative_cell)
  # }
  # print(rownames(representative_cells))
  # first_cells <- data.frame(cluster_id = character(), cell_name = character(), stringsAsFactors = FALSE)
  #
  # for (cluster_id in unique(Idents(obj))) {
  #   cluster_cells <- which(Idents(obj) == cluster_id)
  #   first_cell_index <- cluster_cells[1]
  #   first_cell_name <- colnames(obj)[first_cell_index]
  #
  #   first_cells <- rbind(first_cells, data.frame(cluster_id = cluster_id, cell_name = first_cell_name))
  # }
  #
  # print(first_cells$cell_name)
  # data <- GetAssayData(object = obj, cells = first_cells$cell_name)
  # print(data)
  # distances <- dist(data, method = "euclidean")
  # print(distances)
  #
  # mds_result <- isoMDS(distances, k = 2)
  #
  # plot(mds_result$points, type = "n", xlab = "MDS Dimension 1", ylab = "MDS Dimension 2", main = "MDS Plot of Seurat Clusters")
  # points(mds_result$points, col = rainbow(length(unique(seurat_obj$meta.data$cluster))), pch = 16)
  # legend("topright", legend = unique(seurat_obj$meta.data$cluster), col = rainbow(length(unique(seurat_obj$meta.data$cluster))), pch = 16, title = "Clusters")
}

# visualize annotation
create_annotation_UMAP <- function(obj, col, pc, resolution, values) {
  # annotation list
  annotation <- list(ambiguous = list(0), cd45.1 = list(1, 2), cd45.2 = list(3))

  tryCatch(
  {
    obj <- FindNeighbors(obj, dims = 1:pc)
    obj <- FindClusters(obj, resolution = resolution)

    # create new metadata
    new_metadata <- rep(NA, nrow(obj))

    # assign cluster with annotation
    for (key in names(annotation)) {
      clusters <- unlist(annotation[[key]])
      for (cluster in clusters) {
        new_metadata[obj$seurat_clusters == cluster] <- key
      }
    }

    # add new annotation data to seurat obj
    obj <- AddMetaData(obj, metadata = new_metadata, col.name = "Annotation")
    Run
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
      remove_modal_spinner()
      return(umap)
    }
  )
}

# sankey plot
create_sankey_plot <- function(obj, values, pc, resolution) {
  # obj <- FindNeighbors(obj, dims = 1:pc)
  # obj <- FindClusters(obj, resolution = resolution)
  annotation <- list(ambiguous = list(0), cd45.1 = list(1, 2), cd45.2 = list(3))
  cluster_info <- as.data.frame(obj@active.ident)
  transitions <- data.frame(From = integer(), To = integer())
  # record transition
  for (i in 2:nrow(cluster_info)) {
    from_cluster <- cluster_info[i - 1, 1]
    to_cluster <- cluster_info[i, 1]
    transition <- data.frame(From = from_cluster, To = to_cluster)
    transitions <- rbind(transitions, transition)
  }

  # calculate frequency of transition
  transition_counts <- table(transitions)

  # create transition
  transition_data <- data.frame(transition_counts)
  colnames(transition_data) <- c("From", "To", "Value")
  transition_data$From <- as.integer(as.character(transition_data$From))
  transition_data$To <- as.integer(as.character(transition_data$To))

  # Nodes parameter
  all_nodes <- data.frame(name = c(
    as.character(unique(transition_data$From)),
    as.character(names(annotation))
  ))

  # adjust to list
  num_clusters <- length(unique(obj@active.ident))
  for (i in seq_along(transition_data$To)) {
    for (cluster in 0:num_clusters) {
      if (transition_data$To[i] == cluster) {
        transition_data$To[i] <- num_clusters - 1 + as.numeric(
          which(sapply(annotation, function(x) cluster %in% x))[1]
        )
        break
      }
    }
  }
  transition_data$To <- as.integer(transition_data$To)

  # create sankey
  sankey <- sankeyNetwork(
    Links = transition_data,
    Source = "From",
    Nodes = all_nodes,
    Target = "To",
    Value = "Value",
    units = "Cell Counts",
    width = 400,
    height = 200
  )
  values$sankey <- sankey
  return(sankey)
}
