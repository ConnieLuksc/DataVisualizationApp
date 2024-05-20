source("functions/ClusteringFunc.R")
source("functions/InitialStepFunc.R")

# Define server logic
server <- function(input, output, session) {
  pcValue <- reactiveVal(NA)
  resolutionValue <- reactiveVal(NA)
  maximum_cate_value <- reactiveVal(NA)
  options(shiny.maxRequestSize = 800 * 1024^2)
  clicked_link <- reactiveVal(FALSE)
  values <- reactiveValues()
  values$saved_list <- list()
  values$count <- 2
  values$selected_gene <- list(genes_1 = NULL, genes_2 = NULL, genes_3 = NULL, genes_4 = NULL)
  values$featurePlots <- list(featurePlot_1 = NULL, featurePlot_2 = NULL, featurePlot_3 = NULL, featurePlot_4 = NULL)
  values$violinPlotGenes <- list(violinPlotGene_1 = NULL, violinPlotGene_2 = NULL, violinPlotGene_3 = NULL, violinPlotGene_4 = NULL)
  values$genes <- list(gene_1 = NULL, gene_2 = NULL, gene_3 = NULL, gene_4 = NULL)
  values$gene_1 <- NULL
  values$gene_2 <- NULL
  values$gene_3 <- NULL
  values$gene_4 <- NULL
  values$annotations <- list()
  values$cluster_num <- 0
  values$previous_selected_layer <- ""
  values$annotation_show <- reactive({
    rep("NA", times = values$cluster_num)
  })
  values$layerJudge <- FALSE

  data <- reactiveVal()
  observe({
    req(input$file)
    data(readRDS(input$file$datapath))
    updateSelectInput(session, "layerSelect",
      choices = names(data()@assays$RNA@layers)
    )
  })


  updateFilter <- function(enable = TRUE) {
    if (enable) {
      # Load the Seurat object from the RDS file
      tmp <- load_seurat_obj(input$file$datapath)
      data(tmp)
    } else {
      values$obj <- NULL
    }
  }

  observeEvent(input$layerSelect, {
    req(input$file, input$layerSelect)
    tmp <- data()

    if (input$layerSelect != "" && input$layerSelect != values$previous_selected_layer) {
      # Extract the specific counts layer
      values$previous_selected_layer <- input$layerSelect
      counts_layer <- tmp@assays[["RNA"]]@layers[[input$layerSelect]]

      # Set min.cells as total number of cells in the sample * 0.01
      min.cells <- (ncol(counts_layer) * 0.01) %>% round()

      # Create a new Seurat object using the extracted layer
      obj <- CreateSeuratObject(
        counts = counts_layer,
        min.cells = min.cells,
        min.features = 200
      )
      obj$group <- "all"
      # obj@assays[["RNA"]]@layers[["counts.Gene Expression.Day0_BMC_STIA_Macs"]] <- NULL
      # obj@assays[["RNA"]]@layers[["counts.Gene Expression.Day7_BMC_STIA_Macs"]] <- NULL
      obj <- PercentageFeatureSet(obj, pattern = "^MT-|^MT.|^mt-", col.name = "percent.mt")
      values$filter_violinPlot <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
      output$filter_violinPlot <- renderPlot(values$filter_violinPlot)
      values$feature_scatter <- create_feature_scatter(obj)
      output$feature_scatter <- renderPlot(values$feature_scatter)
      updateNumericInput(session, "feature_upper", value = max(obj@meta.data[["nFeature_RNA"]]))
      updateNumericInput(session, "feature_lower", value = min(obj@meta.data[["nFeature_RNA"]]))
      updateNumericInput(session, "count_upper", value = max(obj@meta.data[["nCount_RNA"]]))
      updateNumericInput(session, "count_lower", value = min(obj@meta.data[["nCount_RNA"]]))
      updateNumericInput(session, "percent_upper", value = max(obj@meta.data[["percent.mt"]]))
      updateNumericInput(session, "percent_lower", value = min(obj@meta.data[["percent.mt"]]))
      values$obj <- obj
    }
  })

  observe({
    updateUI(!is.null(input$file))
    updateFilter(!is.null(input$file))
  })

  # Filter
  observeEvent(input$filter, {
    tryCatch(
      {
        obj <- subset(values$obj, subset = nFeature_RNA > input$feature_lower &
          nFeature_RNA < input$feature_upper &
          nCount_RNA > input$count_lower &
          nCount_RNA < input$count_upper &
          percent.mt > input$percent_lower &
          percent.mt < input$percent_upper)
        values$filter_violinPlot <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "group", ncol = 3, pt.size = 0)
        output$filter_violinPlot <- renderPlot(values$filter_violinPlot)
        values$feature_scatter <- create_feature_scatter(obj)
        output$feature_scatter <- renderPlot(values$feature_scatter)
        values$obj <- obj
      },
      error = function(e) {
        # Error handling
        showModal(modalDialog(
          title = "Error",
          paste("An error has occurred:", e$message),
          easyClose = TRUE,
          footer = NULL
        ))
      }
    )
  })

  # Reset Filter
  observeEvent(input$filter_reset, {
    updateFilter(TRUE)
  })

  observeEvent(input$normalize, {
    # browser()
    tryCatch(
      {
        withProgress(message = "Normalization in progress...", value = 0, {
          obj <- values$obj
          obj <- normalizeData(obj, input$normalization_method, input$using_log, input$norm_parameter, input$sct_parameter)
          obj <- runPCA(obj, input$num_pcs)
          variableFeatures <- findVariableFeatures(obj, input$num_features)
          obj <- variableFeatures$object
          top10 <- variableFeatures$top10

          values$obj <- obj

          values$elbowPlot <- {
            ElbowPlot(obj, ndims = input$num_pcs)
          }
          output$elbowPlot <- renderPlot(values$elbowPlot)
          values$feature_selection <- {
            plotVariableFeatures(obj, top10)
          }
          output$feature_selection <- renderPlot(values$feature_selection)

          for (i in 1:10) {
            Sys.sleep(0.5)
            incProgress(1 / 10, detail = "Normalization complete! Please go to Clustering tab to visualize data")
          }
        })
      },
      error = function(e) {
        showModal(modalDialog(
          title = "Error",
          paste("An error has occurred:", e$message),
          easyClose = TRUE,
          footer = NULL
        ))
      }
    )
  })


  observeEvent(input$run, {
    shinyjs::disable("run")
    shinyjs::disable("update")
    pcValue(as.numeric(input$pc))
    resolutionValue(as.numeric(input$resolution))
    values$annotations <- input$annotation_column
    # Assuming load_seurat_obj is a function you've defined to load the Seurat object
    obj <- values$obj
    values$run_triggered <- reactiveVal(FALSE)

    if (is.vector(values$obj)) {
      # Handle error in file upload or object loading
      showModal(modalDialog(
        title = "Error with file",
        "There is an error with the file you uploaded."
      ))
      shinyjs::enable("run")
    } else {
      if (!is.na(pcValue()) && !is.na(resolutionValue())) {
        show_modal_spinner(text = "Preparing plots...")
        run_umap(obj, pcValue(), resolutionValue(), values)
      }

      output$umap <- renderPlot({
        if (!is.na(pcValue()) && !is.na(resolutionValue()) && !is.null(values$obj)) {
          create_metadata_UMAP(values$obj, values)
        } else {
          ggplot() +
            theme_void() +
            geom_text(aes(x = 0.5, y = 0.5, label = str_wrap("Please input the parameters", width = 20)), size = 12, color = "gray73", fontface = "bold") + # nolint
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
        }
      })

      output$violinPlot <- renderPlot({
        if (!is.na(pcValue()) &&
          !is.na(resolutionValue()) &&
          !is.null(values$obj)) {
          values$obj[["percent.mt"]] <- PercentageFeatureSet(values$obj, pattern = "^MT-")
          violinPlot <- VlnPlot(values$obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
          values$violinPlot <- violinPlot
          violinPlot
        }
      })

      output$normalization <- renderText({
        paste("Normalization Method: ", input$normalization_method)
      })

      output$normalization_parameter <- renderText({
        paste("Normalization Parameter: ", input$parameter)
      })

      output$PC_value <- renderText({
        paste("PC value: ", input$num_pcs)
      })

      output$Number_of_Variable_Features <- renderText({
        paste("Number of Variable Features: ", input$num_features)
      })

      output$heatmapPlot <- renderPlot(
        {
          create_heatmap(values)
        },
        height = function() {
          350
        },
        width = function() {
          500
        }
      )


      output$cluster_cell_counts <- DT::renderDataTable(
        {
          cluster_ids <- Idents(values$obj)
          values$cluster_cell_counts <- table(cluster_ids)
          values$cluster_num <- length(names(values$cluster_cell_counts))
          data.frame(
            Cluster = names(values$cluster_cell_counts),
            Count = as.numeric(values$cluster_cell_counts)
          )
        },
        options = list(paging = FALSE),
        rownames = FALSE
      )

      output$mdsPlot <- renderPlot({
        if (!is.null(values$obj)) {
          create_mds_plot(values$obj, values)
        }
      })

      output$umap_annotation <- renderPlot({
        create_annotation_UMAP(values$obj, pcValue(), resolutionValue(), values, values$annotations)
      })

      output$unavailable_sankey <- renderPlot({
        ggplot() +
          theme_void() +
          geom_text(
            aes(
              x = 0.5, y = 0.5,
              label = str_wrap("Annotations missing", width = 20)
            ),
            size = 12, color = "gray73", fontface = "bold"
          ) +
          theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
      })

      output$sankeyPlot <- renderUI({
        create_sankey_plot(values$obj, values, pcValue(), resolutionValue(), values$annotations, values$cluster_cell_counts)
      })


      output$text_output <- renderText({
        remove_modal_spinner()
        return("Current #Cells/Cluster")
      })

      shinyjs::enable("run")
      shinyjs::enable("save")
      shinyjs::enable("update")
    }
  })


  # Clear all sidebar inputs when 'Reset' button is clicked
  observeEvent(input$reset, {
    # Reset file input
    shinyjs::reset("file")

    # Clear plots and disable the 'Run' button
    output$umap <- renderPlot({
      NULL
    })
    output$featurePlot <- renderPlot({
      NULL
    })
    output$violinPlot <- renderPlot({
      NULL
    })
    output$violinPlotGene <- renderPlot({
      NULL
    })
    output$mdsPlot <- renderPlot({
      NULL
    })

    shinyjs::disable("run")
  })


  # monitor change of saved_list
  observeEvent(values$saved_list, {
    if (!is.null(values$saved_list)) {
      current_saved_list <- values$saved_list

      # update saved area part
      output$dynamic_elements <- renderUI({
        dynamic_elements <- lapply(
          seq_along(current_saved_list),
          function(key) {
            fluidRow(
              column(
                width = 12,
                actionButton(
                  inputId = paste0("button_", key),
                  label = paste(
                    "PC: ", current_saved_list[[key]]$pc,
                    " Resolution: ", current_saved_list[[key]]$resolution,
                    "Gene: ", current_saved_list[[key]]$gene
                  )
                ),
                DT::dataTableOutput(outputId = paste0("verbatim_output_", key))
              )
            )
          }
        )
        # add margin
        div_with_margin <- lapply(dynamic_elements, function(elem) {
          div(style = "margin-top: 20px;", elem)
        })

        do.call(tagList, div_with_margin)
        # do.call(tagList, dynamic_elements) # nolint
      })

      lapply(seq_along(current_saved_list), function(key) {
        output[[paste0("verbatim_output_", key)]] <- DT::renderDataTable(
          {
            data.frame(
              Cluster = names(current_saved_list[[key]]$cluster),
              Count = as.numeric(current_saved_list[[key]]$cluster)
            )
          },
          options = list(paging = FALSE),
          rownames = FALSE
        )

        observeEvent(input[[paste0("button_", key)]],
          {
            loaded_seurat <- load_seurat_obj(current_saved_list[[key]]$file)
            if (!is.vector(loaded_seurat)) {
              # values$obj <- loaded_seurat # nolint
              output[[paste0("pc1", values$count)]] <- renderText({
                paste("PC: ", current_saved_list[[key]]$pc)
              })

              output[[paste0("resolution1", values$count)]] <- renderText({
                paste("Resolution: ", current_saved_list[[key]]$resolution)
              })
              values$gene_1 <- current_saved_list[[key]]$gene_1
              output[[paste0("gene1", values$count)]] <- renderText({
                paste("Gene1: ", values$gene_1)
              })
              values$gene_2 <- current_saved_list[[key]]$gene_2
              output[[paste0("gene2", values$count)]] <- renderText({
                paste("Gene2: ", values$gene_2)
              })
              values$gene_3 <- current_saved_list[[key]]$gene_3
              output[[paste0("gene3", values$count)]] <- renderText({
                paste("Gene3: ", values$gene_3)
              })
              values$gene_4 <- current_saved_list[[key]]$gene_4
              output[[paste0("gene4", values$count)]] <- renderText({
                paste("Gene4: ", values$gene_4)
              })
              values$annotations <- current_saved_list[[key]]$annotation
              output[[paste0("annotation1", values$count)]] <- renderText({
                paste("Annotation: ", values$annotations)
              })
              output[[paste0("umap", values$count)]] <- renderPlot(current_saved_list[[key]]$umap) # nolint
              output[[paste0("violinPlot", values$count)]] <- renderPlot(current_saved_list[[key]]$violin) # nolint
              output[[paste0("featurePlot_1", values$count)]] <- renderPlot(current_saved_list[[key]]$featurePlot_1) # nolint
              output[[paste0("violinPlotGene_1", values$count)]] <- renderPlot(current_saved_list[[key]]$violinPlotGene_1) # nolint
              output[[paste0("featurePlot_2", values$count)]] <- renderPlot(current_saved_list[[key]]$featurePlot_2) # nolint
              output[[paste0("violinPlotGene_2", values$count)]] <- renderPlot(current_saved_list[[key]]$violinPlotGene_2) # nolint
              output[[paste0("featurePlot_3", values$count)]] <- renderPlot(current_saved_list[[key]]$featurePlot_3) # nolint
              output[[paste0("violinPlotGene_3", values$count)]] <- renderPlot(current_saved_list[[key]]$violinPlotGene_3) # nolint
              output[[paste0("featurePlot_4", values$count)]] <- renderPlot(current_saved_list[[key]]$featurePlot_4) # nolint
              output[[paste0("violinPlotGene_4", values$count)]] <- renderPlot(current_saved_list[[key]]$violinPlotGene_4) # nolint
              output[[paste0("sankeyPlot", values$count)]] <- renderUI(current_saved_list[[key]]$sankey) # nolint
              output[[paste0("heatmapPlot", values$count)]] <- renderPlot(current_saved_list[[key]]$heatmap)
              output[[paste0("mdsPlot", values$count)]] <- renderPlot({
                plotMDS(current_saved_list[[key]]$plotmds, col = current_saved_list[[key]]$mds_color, pch = 20, cex = 2)
                title("MDS Plot")
                legend("topright", legend = colnames(current_saved_list[[key]]$mds_legend), col = current_saved_list[[key]]$mds_color, pch = 20, cex = 0.8, pt.cex = 0.8, title = "Group")
              })
              output[[paste0("normalization", values$count)]] <- renderText({
                paste("Normalization Method: ", input$normalization_method)
              })
              output[[paste0("normalization_parameter", values$count)]] <- renderText({
                paste("Normalization Parameter: ", input$parameter)
              })
              output[[paste0("PC_value", values$count)]] <- renderText({
                paste("PC value: ", input$num_pcs)
              })
              output[[paste0("Number_of_Variable_Features", values$count)]] <- renderText({
                paste("Number of Variable Features: ", input$num_features)
              })
              output[[paste0("umap_annotation", values$count)]] <- renderPlot(current_saved_list[[key]]$annotation_umap) # nolint
              output[[paste0("textoutput", values$count)]] <- renderText({
                return("Current #Cells/Cluster")
              })
              output[[paste0("cluster_cell_counts", values$count)]] <- DT::renderDataTable( # nolint
                {
                  data.frame(
                    Cluster = names(current_saved_list[[key]]$cluster),
                    Count = as.numeric(current_saved_list[[key]]$cluster)
                  )
                },
                options = list(paging = FALSE),
                rownames = FALSE
              )
              insertTab(
                inputId = "main_tabs",
                tabPanel(
                  paste("Tab", values$count),
                  fluidRow(
                    column(
                      2,
                      verbatimTextOutput(paste0("pc1", values$count)),
                      verbatimTextOutput(paste0("resolution1", values$count)),
                      verbatimTextOutput(paste0("gene1", values$count)),
                      verbatimTextOutput(paste0("gene2", values$count)),
                      verbatimTextOutput(paste0("gene3", values$count)),
                      verbatimTextOutput(paste0("gene4", values$count)),
                      verbatimTextOutput(paste0("annotation1", values$count)),
                      verbatimTextOutput(paste0("normalization", values$count)),
                      verbatimTextOutput(paste0("normalization_parameter", values$count)),
                      verbatimTextOutput(paste0("PC_value", values$count)),
                      verbatimTextOutput(paste0("Number_of_Variable_Features", values$count)),
                      downloadButton("download_r", "Download R object")
                    ),
                    column(
                      8,
                      fluidRow(
                        column(
                          5,
                          plotOutput(paste0("violinPlot", values$count))
                        ),
                        column(
                          7,
                          plotOutput(paste0("umap", values$count))
                        ),
                      ),
                      fluidRow(
                        column(3, plotOutput(paste0("featurePlot_1", values$count), height = "225px")),
                        column(3, plotOutput(paste0("featurePlot_2", values$count), height = "225px")),
                        column(3, plotOutput(paste0("featurePlot_3", values$count), height = "225px")),
                        column(3, plotOutput(paste0("featurePlot_4", values$count), height = "225px"))
                      ),
                      fluidRow(
                        column(3, plotOutput(paste0("violinPlotGene_1", values$count), height = "225px")),
                        column(3, plotOutput(paste0("violinPlotGene_2", values$count), height = "225px")),
                        column(3, plotOutput(paste0("violinPlotGene_3", values$count), height = "225px")),
                        column(3, plotOutput(paste0("violinPlotGene_4", values$count), height = "225px"))
                      ),
                      fluidRow(
                        column(
                          6,
                          plotOutput(paste0("heatmapPlot", values$count))
                        ),
                        column(
                          6,
                          plotOutput(paste0("mdsPlot", values$count))
                        )
                      ),
                      fluidRow(
                        column(
                          6,
                          plotOutput(paste0("umap_annotation", values$count))
                        ),
                        column(
                          6,
                          uiOutput(paste0("sankeyPlot", values$count))
                        )
                      )
                    ),
                    column(
                      width = 2,
                      style = "overflow-y: scroll; max-height: 600px;",
                      textOutput(paste0("textoutput", values$count)),
                      DT::dataTableOutput(
                        paste0("cluster_cell_counts", values$count)
                      ),
                    )
                  )
                ),
                select = TRUE,
              )
              shinyjs::runjs('$("#pc1").prop("disabled", true);')
              shinyjs::runjs('$("#gene1").prop("disabled", true);')
              shinyjs::runjs('$("#resolution1").prop("disabled", true);')
              values$count <- values$count + 1
            }
          },
          ignoreInit = TRUE
        )
      })
    }
  })

  # update pc value
  observeEvent(input$update, {
    maximum_cate_value(as.numeric(input$maximum_categories))
    values$annotations <- input$annotation_column
    values$genes_1 <- input$genes_1
    values$gene_1 <- input$genes_1
    values$genes_2 <- input$genes_2
    values$gene_2 <- input$genes_2
    values$genes_3 <- input$genes_3
    values$gene_3 <- input$genes_3
    values$genes_4 <- input$genes_4
    values$gene_4 <- input$genes_4

    if (pcValue() != input$pc || resolutionValue() != input$resolution) {
      output$featurePlot_1 <- renderPlot(create_feature_plot(values$obj, values$genes_1, values))
      output$featurePlot_2 <- renderPlot(create_feature_plot(values$obj, values$genes_2, values))
      output$featurePlot_3 <- renderPlot(create_feature_plot(values$obj, values$genes_3, values))
      output$featurePlot_4 <- renderPlot(create_feature_plot(values$obj, values$genes_4, values))
      output$violinPlotGene_1 <- renderPlot(create_violin_plot(values$obj, values$genes_1, values, ncol = NULL, pt.size = 0))
      output$violinPlotGene_2 <- renderPlot(create_violin_plot(values$obj, values$genes_2, values, ncol = NULL, pt.size = 0))
      output$violinPlotGene_3 <- renderPlot(create_violin_plot(values$obj, values$genes_3, values, ncol = NULL, pt.size = 0))
      output$violinPlotGene_4 <- renderPlot(create_violin_plot(values$obj, values$genes_4, values, ncol = NULL, pt.size = 0))
    }
    pcValue(as.numeric(input$pc))
    resolutionValue(as.numeric(input$resolution))

    create_plot_output <- function(i, gene_index, umapped_obj, gene) {
      output[[paste0("featurePlot_", gene_index, i)]] <- renderPlot({
        create_feature_plot(umapped_obj, gene, values)
      })
      output[[paste0("violinPlotGene_", gene_index, i)]] <- renderPlot({
        create_violin_plot(umapped_obj, gene, values, ncol = NULL, pt.size = 0)
      })
    }

    if (length(values$saved_list) > 0) {
      for (i in seq(from = 2, to = values$count - 1)) {
        local({
          current_saved_list <- values$saved_list
          umapped_obj <- current_saved_list[[i - 1]]$umapped_obj
          pc <- current_saved_list[[i - 1]]$pc
          resolution <- current_saved_list[[i - 1]]$resolution
          create_plot_output(i, 1, umapped_obj, values$gene_1)
          create_plot_output(i, 2, umapped_obj, values$gene_2)
          create_plot_output(i, 3, umapped_obj, values$gene_3)
          create_plot_output(i, 4, umapped_obj, values$gene_4)
          output[[paste0("umap_annotation", i)]] <- renderPlot({
            create_annotation_UMAP(umapped_obj, pc, resolution, values, values$annotations)
          })
          output[[paste0("sankeyPlot", i)]] <- renderUI({
            create_sankey_plot(umapped_obj, values, pc, resolution, values$annotations, values$cluster_cell_counts)
          })
        })
      }
    }
  })

  observeEvent(maximum_cate_value(), {
    if (!is.null(values$obj)) {
      max_categories <- input$maximum_categories
      all_columns <- colnames(values$obj@meta.data)

      selected_columns <- c()

      for (column in all_columns) {
        if (length(unique(values$obj@meta.data[[column]])) <= max_categories) {
          selected_columns <- c(selected_columns, column)
        }
      }
      updateSelectizeInput(session, "annotation_column", choices = selected_columns)
    }
  })

  observeEvent(input$save, {
    new_index <- paste0("item", length(values$saved_list) + 1)
    saved_list_tmp <- list(
      file = input$file$datapath, pc = pcValue(), umapped_obj = values$umapped_obj,
      resolution = resolutionValue(), gene = input$gene, annotation = values$annotations,
      cluster = values$cluster_cell_counts, umap = values$umap,
      violin = values$violinPlot, heatmap = values$heatmap,
      featurePlot_1 = values$featurePlots[["featurePlot_1"]], violinPlotGene_1 = values$violinPlotGenes[["violinPlotGene_1"]],
      featurePlot_2 = values$featurePlots[["featurePlot_2"]], violinPlotGene_2 = values$violinPlotGenes[["violinPlotGene_2"]],
      featurePlot_3 = values$featurePlots[["featurePlot_3"]], violinPlotGene_3 = values$violinPlotGenes[["violinPlotGene_3"]],
      featurePlot_4 = values$featurePlots[["featurePlot_4"]], violinPlotGene_4 = values$violinPlotGenes[["violinPlotGene_4"]],
      geneViolin = values$violin, annotation_umap = values$umap_annotation,
      sankey = values$sankey, plotmds = values$mds, mds_legend = values$mds_legend, mds_color = values$mds_color,
      gene_1 = values$genes[["gene_1"]], gene_2 = values$genes[["gene_2"]],
      gene_3 = values$genes[["gene_3"]], gene_4 = values$genes[["gene_4"]]
    )
    values$saved_list[[new_index]] <- saved_list_tmp
  })

  plots <- reactiveValues()
  output$download_initial <- downloadHandler(
    filename = function() {
      paste("plots_initial_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, onefile = TRUE, width = 15, height = 9)
      line1 <- sprintf("nFeature RNA: %d - %d", input$feature_lower, input$feature_upper)
      line2 <- sprintf("nCount RNA: %d - %d", input$count_lower, input$count_upper)
      line3 <- sprintf("Percent mt: %.2f - %.2f", input$percent_lower, input$percent_upper)
      plots$filter_violinPlot <- values$filter_violinPlot +
        theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))
      print(plots$filter_violinPlot)
      grid.text(line1, x = 0.5, y = 0.05, just = "bottom", gp = gpar(fontsize = 15))
      grid.text(line2, x = 0.5, y = 0.03, just = "bottom", gp = gpar(fontsize = 15))
      grid.text(line3, x = 0.5, y = 0.01, just = "bottom", gp = gpar(fontsize = 15))

      plots$feature_scatter <- values$feature_scatter +
        theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))
      print(plots$feature_scatter)
      grid.text(line1, x = 0.5, y = 0.05, just = "bottom", gp = gpar(fontsize = 15))
      grid.text(line2, x = 0.5, y = 0.03, just = "bottom", gp = gpar(fontsize = 15))
      grid.text(line3, x = 0.5, y = 0.01, just = "bottom", gp = gpar(fontsize = 15))

      plots$elbowPlot <- values$elbowPlot +
        theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))
      print(plots$elbowPlot)
      grid.text(paste("PC value: ", input$num_pcs), x = 0.5, y = 0.05, just = "bottom", gp = gpar(fontsize = 15))

      plots$feature_selection <- values$feature_selection +
        theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))
      print(plots$feature_selection)
      grid.text(paste("Number of Variable Features:", input$num_features), x = 0.5, y = 0.05, just = "bottom", gp = gpar(fontsize = 15))

      dev.off()
    }
  )


  output$download_clustering <- downloadHandler(
    filename = function() {
      paste("plots_clustering_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, onefile = TRUE, width = 15, height = 9)
      plots$violinPlot <- values$violinPlot +
        theme(
          plot.margin = unit(c(2, 2, 2, 2), "cm"),
          plot.title = element_text(hjust = 0.5))
      print(plots$violinPlot)
      grid.text(paste("PC: ", input$pc), x = 0.5, y = 0.05, just = "bottom", gp = gpar(fontsize = 15))
      grid.text(paste("Resolution: ", input$resolution), x = 0.5, y = 0.03, just = "bottom", gp = gpar(fontsize = 15))

      plots$umap <- values$umap +
        ggtitle("UMAP") +
        theme(
          plot.margin = unit(c(2, 2, 2, 2), "cm"),
          plot.title = element_text(hjust = 0.5))
      print(plots$umap)
      grid.text(paste("PC: ", input$pc), x = 0.5, y = 0.05, just = "bottom", gp = gpar(fontsize = 15))
      grid.text(paste("Resolution: ", input$resolution), x = 0.5, y = 0.03, just = "bottom", gp = gpar(fontsize = 15))

      plots$umap_annotation <- values$umap_annotation + 
        ggtitle("UMAP Annotations") + theme(plot.title = element_text(hjust = 0.5))
      print(plots$umap_annotation)

      plotMDS(values$mds, col = values$mds_color, pch = 20, cex = 2)
      title("MDS Plot")
      legend("topright", legend = colnames(values$mds_legend), col = values$mds_color, pch = 20, cex = 0.8, pt.cex = 0.8, title = "Group")

      grid.newpage()
      print(values$heatmap)
      print(values$featurePlots[["featurePlot_1"]])
      print(values$violinPlotGenes[["violinPlotGene_1"]])
      print(values$featurePlots[["featurePlot_2"]])
      print(values$violinPlotGenes[["violinPlotGene_2"]])
      print(values$featurePlots[["featurePlot_3"]])
      print(values$violinPlotGenes[["violinPlotGene_3"]])
      print(values$featurePlots[["featurePlot_4"]])
      print(values$violinPlotGenes[["violinPlotGene_4"]])
      dev.off()
    }
  )

  output$download_r <- downloadHandler(
    filename = function() {
      "seurat_obj.rds"
    },
    content = function(file) {
      saveRDS(values$obj, file)
    }
  )


  # choose annotation colum from meta data
  # observeEvent(input$annotation_column, {
  #   # print(input$annotation_column)
  #   values$annotations <- input$annotation_column
  # })


  # feature plot and violin plot for gene
  observeEvent(values$genes_1, {
    values$featurePlots[["featurePlot_1"]] <- create_feature_plot(values$obj, values$genes_1, values)
    output$featurePlot_1 <- renderPlot(values$featurePlots[["featurePlot_1"]])
    values$violinPlotGenes[["violinPlotGene_1"]] <- create_violin_plot(values$obj, values$genes_1, values, ncol = NULL, pt.size = 0)
    values$genes[["gene_1"]] <- values$genes_1
    output$violinPlotGene_1 <- renderPlot(values$violinPlotGenes[["violinPlotGene_1"]])
    values$selected_genes[["genes_1"]] <- values$genes_1
  })
  observeEvent(values$genes_2, {
    values$featurePlots[["featurePlot_2"]] <- create_feature_plot(values$obj, values$genes_2, values)
    output$featurePlot_2 <- renderPlot(values$featurePlots[["featurePlot_2"]])
    values$violinPlotGenes[["violinPlotGene_2"]] <- create_violin_plot(values$obj, values$genes_2, values, ncol = NULL, pt.size = 0)
    values$genes[["gene_2"]] <- values$genes_2
    output$violinPlotGene_2 <- renderPlot(values$violinPlotGenes[["violinPlotGene_2"]])
    values$selected_genes[["genes_2"]] <- values$genes_2
  })
  observeEvent(values$genes_3, {
    values$featurePlots[["featurePlot_3"]] <- create_feature_plot(values$obj, values$genes_3, values)
    output$featurePlot_3 <- renderPlot(values$featurePlots[["featurePlot_3"]])
    values$violinPlotGenes[["violinPlotGene_3"]] <- create_violin_plot(values$obj, values$genes_3, values, ncol = NULL, pt.size = 0)
    values$genes[["gene_3"]] <- values$genes_3
    output$violinPlotGene_3 <- renderPlot(values$violinPlotGenes[["violinPlotGene_3"]])
    values$selected_genes[["genes_3"]] <- values$genes_3
  })
  observeEvent(values$genes_4, {
    values$featurePlots[["featurePlot_4"]] <- create_feature_plot(values$obj, values$genes_4, values)
    output$featurePlot_4 <- renderPlot(values$featurePlots[["featurePlot_4"]])
    values$violinPlotGenes[["violinPlotGene_4"]] <- create_violin_plot(values$obj, values$genes_4, values, ncol = NULL, pt.size = 0)
    values$genes[["gene_4"]] <- values$genes_4
    output$violinPlotGene_4 <- renderPlot(values$violinPlotGenes[["violinPlotGene_4"]])
    values$selected_genes[["genes_4"]] <- values$genes_4
  })

  observeEvent(values$obj, {
    if (!is.null(values$obj)) {
      updateSelectizeInput(session, "genes_1", choices = rownames(values$obj), selected = values$selected_genes[["genes_1"]], server = TRUE)
      updateSelectizeInput(session, "genes_2", choices = rownames(values$obj), selected = values$selected_genes[["genes_2"]], server = TRUE)
      updateSelectizeInput(session, "genes_3", choices = rownames(values$obj), selected = values$selected_genes[["genes_3"]], server = TRUE)
      updateSelectizeInput(session, "genes_4", choices = rownames(values$obj), selected = values$selected_genes[["genes_4"]], server = TRUE)
      # updateSelectizeInput(session, "annotation_column", choices = colnames(values$obj@meta.data))
      all_columns <- colnames(values$obj@meta.data)
      selected_columns <- c()
      for (column in all_columns) {
        if (length(unique(values$obj@meta.data[[column]])) <= 100) {
          selected_columns <- c(selected_columns, column)
        }
      }
      updateSelectizeInput(session, "annotation_column", choices = selected_columns)
      updateSelectInput(session, "clusters", choices = names(values$cluster_cell_counts))
    }
  })
}
