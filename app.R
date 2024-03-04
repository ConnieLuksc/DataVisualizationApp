source('global.R')

# Define UI
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Single-Cell Clustering"),
  tabsetPanel(
    id = "main_tabs",
    tabPanel(
      "Tab 1",
      fluidRow(
        column(2,
               fileInput("file", "Upload File", multiple = TRUE, accept = c('.rds')),
               actionButton("reset", "Reset", class = "reset-btn"),
               actionButton("run", "Run", class = "run-btn"),
               numericInput("pc", "PC", value = NA),
               numericInput("resolution", "Resolution", value = NA, step = 0.1),
               selectizeInput('genes_1', 'Genes', choices = NULL, multiple = TRUE),
               selectizeInput('genes_2', NULL, choices = NULL, multiple = TRUE),
               selectizeInput('genes_3', NULL, choices = NULL, multiple = TRUE),
               selectizeInput('genes_4', NULL, choices = NULL, multiple = TRUE),
               selectizeInput('annotation_column', label = 'annotation_column',  choices = NULL),
               selectInput('clusters', 'Clusters', choices = NULL, multiple = TRUE, selectize = TRUE),
               textInput("annotation", "Annotation"),
               selectizeInput("gene", "Genes", choices = NULL),
               actionButton("annotate", "Annotate"),
               actionButton("save", "Save", class = "save-btn")
        ),
        column(
          8,
          fluidRow(
            column(5, plotOutput("violinPlot")),
            column(7, plotOutput("umap")),
          ),
          fluidRow(
            column(3, plotOutput(outputId = 'featurePlot_1', height = '225px')),
            column(3, plotOutput(outputId = 'featurePlot_2', height = '225px')),
            column(3, plotOutput(outputId = 'featurePlot_3', height = '225px')),
            column(3, plotOutput(outputId = 'featurePlot_4', height = '225px'))
          ),
          fluidRow(
            column(3, plotOutput(outputId = 'violinPlotGene_1', height = '225px')),
            column(3, plotOutput(outputId = 'violinPlotGene_2', height = '225px')),
            column(3, plotOutput(outputId = 'violinPlotGene_3', height = '225px')),
            column(3, plotOutput(outputId = 'violinPlotGene_4', height = '225px'))
          ),
          fluidRow(
            column(
              6,
              plotOutput("heatmapPlot")
            ),
            column(
              6,
              plotOutput(outputId = "mdsPlot")
            )

          ),
          fluidRow(
            column(
              6,
              plotOutput(outputId = "umap_annotation")
            ),
            column(
              6,
              uiOutput(outputId = "sankeyPlot")
            )
          )
        ),
        column(
          width = 2,
          style = "overflow-y: scroll; max-height: 600px;",
          textOutput(outputId = "text_output"),
          # verbatimTextOutput("cluster_cell_counts"),
          DT::dataTableOutput("cluster_cell_counts"),
          uiOutput("dynamic_elements")
        )
      ),
    )
  )
)


# Define server logic
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 300 * 1024^2)
  clicked_link <- reactiveVal(FALSE)
  values <- reactiveValues()
  values$saved_list <- list()
  ignore_button_clicked <- FALSE
  values$count <- 2
  values$selected_gene <- list(genes_1 = NULL, genes_2 = NULL, genes_3 = NULL, genes_4 = NULL)
  values$featurePlots <- list(featurePlot_1 = NULL, featurePlot_2 = NULL, featurePlot_3 = NULL, featurePlot_4 = NULL)
  values$violinPlotGenes <- list(violinPlotGene_1 = NULL, violinPlotGene_2 = NULL, violinPlotGene_3 = NULL, violinPlotGene_4 = NULL)
  values$annotations <- list()
  values$cluster_num <- 0
  values$annotation_show <- reactive({ rep("NA", times = values$cluster_num) })


  updateUI <- function(enable = TRUE) {
    if (enable) {
      shinyjs::enable("run")
      shinyjs::enable("pc")
      shinyjs::enable("resolution")
    } else {
      shinyjs::disable("run")
      shinyjs::disable("pc")
      shinyjs::disable("resolution")
      shinyjs::disable("save")
    }
  }

  observe({
    updateUI(!is.null(input$file))
  })

  observeEvent(input$run, {
    shinyjs::disable("run")

    # Assuming load_seurat_obj is a function you've defined to load the Seurat object
    obj <- load_seurat_obj(input$file$datapath)
    values$obj <- obj
    values$run_triggered <- reactiveVal(FALSE)
    # values$selected_gene <- rownames(obj)[1]

    if (is.vector(values$obj)) {
      # Handle error in file upload or object loading
      showModal(modalDialog(
        title = "Error with file",
        "There is an error with the file you uploaded."
      ))
      shinyjs::enable("run")
    } else {
      output$umap <- renderPlot({
        if (!is.na(input$pc) && !is.na(input$resolution)) {
          show_modal_spinner(text = "Preparing plots...")
          create_metadata_UMAP(
            obj, "seurat_clusters",
            input$pc, input$resolution,
            values
          )
        } else {
          ggplot() +
            theme_void() +
            geom_text(aes(x = 0.5, y = 0.5, label = str_wrap("Please input the parameters", width = 20)), size = 12, color = "gray73", fontface = "bold") + # nolint
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
        }
      })

      # output$featurePlot <- renderPlot({
      #   if (!is.null(input$genes)) {
      #     create_feature_plot(values$obj, input$genes, values)
      #   }
      # })

      output$violinPlot <- renderPlot({
<<<<<<< HEAD
        if (!is.na(input$pc) &&
          !is.na(input$resolution) &&
          !is.null(values$obj)) {
          values$obj[["percent.mt"]] <- PercentageFeatureSet(values$obj, pattern = "^MT-")
          violinPlot <- VlnPlot(values$obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
          values$violinPlot <- violinPlot
          violinPlot
        }
=======
          if (!is.na(input$pc) && !is.na(input$resolution) && !is.null(values$obj)) {
              values$obj[["percent.mt"]] <- PercentageFeatureSet(values$obj, pattern = "^MT-")
              create_violin_plot(values$obj, values = values, ncol = 3, pt.size = 0)
          }
>>>>>>> 2a24397393fc56dc75ab34231bbd069f26fffb1b
      })


      output$violinPlotGene <- renderPlot({
        if (!is.null(input$gene)) {
          create_violin_plot(values$obj, input$gene, values,
                             ncol = NULL, pt.size = 0
          )
        }
      })


      output$heatmapPlot <- renderPlot({
        req(values$obj)
        # Obtain variable genes
        variable_genes <- Seurat::VariableFeatures(values$obj)
        if (is.null(variable_genes) || length(variable_genes) == 0) {
          stop("Variable features not found or the list is empty.")
        }
        # Aggregate expression data for variable genes
        avg_expression <- AggregateExpression(values$obj, features = variable_genes, return.seurat = TRUE)
        data_matrix <- GetAssayData(avg_expression, slot = "data")
        if (!is.matrix(data_matrix) || !is.numeric(data_matrix)) {
          stop("The data matrix is not numeric.")
        }
        # Determine clusters based on the dynamic input
        cluster_assignments <- Idents(values$obj)
        # Find all unique cluster IDs
        all_clusters <- unique(cluster_assignments)
        all_clusters <- sort(all_clusters)
        print("all_clusters")
        print(all_clusters)
        colnames(data_matrix) <- gsub("^g", "", colnames(data_matrix))
        print(paste("Unique cluster IDs:", paste(all_clusters, collapse = ", ")))
        print(paste("Column names in data_matrix:", paste(colnames(data_matrix), collapse = ", ")))
        # Check that data_matrix is not empty after potential subsetting (if needed)
        if (ncol(data_matrix) == 0) {
          stop("The data matrix has no columns.")
        }
        # Calculate the standard deviation for each gene and filter out the genes with zero standard deviation
        non_zero_variance_genes <- apply(data_matrix, 1, var, na.rm = TRUE) > 0
        data_matrix <- data_matrix[non_zero_variance_genes,]
        # Check that data_matrix is not empty after filtering for non-zero variance genes
        if (nrow(data_matrix) == 0) {
          stop("No variable genes found with non-zero variance.")
        }

        # cluster_colors <- grDevices::rainbow(length(all_clusters))
        set.seed(123)
        colors <- rainbow(length(all_clusters))
        print("colors")
        print(colors)
        names(colors) <- all_clusters
        cluster_annotation <- data.frame(Cluster = colnames(data_matrix))
        rownames(cluster_annotation) <- colnames(data_matrix)
        print("cluster_annotation")
        print(head(cluster_annotation))
        annotation_colors = list(Cluster = colors)


        # Calculate the correlation matrix on the subsetted data, handling any remaining NAs
        correlation_matrix <- cor(data_matrix, use = "pairwise.complete.obs")

        print("Inspecting first few rows of correlation_matrix:")
        print(head(correlation_matrix))
        # Plot the heatmap

        values$heatmap <- pheatmap(correlation_matrix,
<<<<<<< HEAD
                                   clustering_distance_rows = "euclidean",
                                   clustering_distance_cols = "euclidean",
                                   clustering_method = "complete",
                                   color = colorRampPalette(c("yellow", "orange", "red"))(50),
                                   annotation_col = cluster_annotation,
                                   annotation_colors = annotation_colors)
=======
                          clustering_distance_rows = "euclidean",
                          clustering_distance_cols = "euclidean",
                          clustering_method = "complete",
                          color = colorRampPalette(c("yellow", "orange", "red"))(50),
                          annotation_col = cluster_annotation,
                          annotation_colors = annotation_colors)

>>>>>>> 2a24397393fc56dc75ab34231bbd069f26fffb1b
        values$heatmap
      }, height = function() { 350 }, width = function() { 500 })



      output$cluster_cell_counts <- DT::renderDataTable({
        cluster_ids <- Idents(values$obj)
        values$cluster_cell_counts <- table(cluster_ids)
        values$cluster_num <- length(names(values$cluster_cell_counts))
        data.frame(
          Cluster = names(values$cluster_cell_counts),
          # Annotation = values$annotation_show,
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
        create_annotation_UMAP(obj, "seurat_clusters", input$pc, input$resolution, values, values$annotations)
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
        create_sankey_plot(values$obj, values, input$pc, input$resolution, values$annotations, values$cluster_cell_counts)
      })


      output$text_output <- renderText({
        remove_modal_spinner()
        return("Current #Cells/Cluster")
      })

      shinyjs::enable("run")
      shinyjs::enable("save")
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

        observeEvent(input[[paste0("button_", key)]], {
          print("button")
          print(ignore_button_clicked)
          if (!ignore_button_clicked) {
            print("loaded")
            loaded_seurat <- load_seurat_obj(current_saved_list[[key]]$file)
            if (!is.vector(loaded_seurat)) {
              # values$obj <- loaded_seurat # nolint
              output[[paste0("umap", values$count)]] <- renderPlot(current_saved_list[[key]]$umap) # nolint
              output[[paste0("violinPlot", values$count)]] <- renderPlot(current_saved_list[[key]]$violin) # nolint
              output[[paste0("featurePlot_1", values$count)]] <- current_saved_list[[key]]$featurePlot_1 # nolint
              output[[paste0("violinPlotGene_1", values$count)]] <- current_saved_list[[key]]$violinPlotGene_1 # nolint
              output[[paste0("featurePlot_2", values$count)]] <- current_saved_list[[key]]$featurePlot_2 # nolint
              output[[paste0("violinPlotGene_2", values$count)]] <- current_saved_list[[key]]$violinPlotGene_2 # nolint
              output[[paste0("featurePlot_3", values$count)]] <- current_saved_list[[key]]$featurePlot_3 # nolint
              output[[paste0("violinPlotGene_3", values$count)]] <- current_saved_list[[key]]$violinPlotGene_3 # nolint
              output[[paste0("featurePlot_4", values$count)]] <- current_saved_list[[key]]$featurePlot_4 # nolint
              output[[paste0("violinPlotGene_4", values$count)]] <- current_saved_list[[key]]$violinPlotGene_4 # nolint
              output[[paste0("sankeyPlot", values$count)]] <- renderUI(current_saved_list[[key]]$sankey) # nolint
              output[[paste0("heatmapPlot", values$count)]] <- renderPlot(current_saved_list[[key]]$heatmap)
              output[[paste0("mdsPlot", values$count)]] <- renderPlot(current_saved_list[[key]]$plotmds)
              print(current_saved_list[[key]]$plotmds)
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
                      numericInput("pc1", "PC",
                                   value = current_saved_list[[key]]$pc
                      ),
                      numericInput("resolution1", "Resolution",
                                   value = current_saved_list[[key]]$resolution, step = 0.1
                      ),
                      selectizeInput("gene1", "Genes",
                                     choices = current_saved_list[[key]]$gene
                      )
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
                        column(3, plotOutput(paste0("featurePlot_1", values$count), height = '225px')),
                        column(3, plotOutput(paste0("featurePlot_2", values$count), height = '225px')),
                        column(3, plotOutput(paste0("featurePlot_3", values$count), height = '225px')),
                        column(3, plotOutput(paste0("featurePlot_4", values$count), height = '225px'))
                      ),
                      fluidRow(
                        column(3, plotOutput(paste0("violinPlotGene_1", values$count), height = '225px')),
                        column(3, plotOutput(paste0("violinPlotGene_2", values$count), height = '225px')),
                        column(3, plotOutput(paste0("violinPlotGene_3", values$count), height = '225px')),
                        column(3, plotOutput(paste0("violinPlotGene_4", values$count), height = '225px'))
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
          }
          ignore_button_clicked <<- FALSE
        })
      })
    }
  })

  observeEvent(input$save, {
    ignore_button_clicked <<- TRUE
    print("save")
    print(ignore_button_clicked)
    new_index <- paste0("item", length(values$saved_list) + 1)
    saved_list_tmp <- list(
      file = input$file$datapath, pc = input$pc,
      resolution = input$resolution, gene = input$gene,
      cluster = values$cluster_cell_counts, umap = values$umap,
      violin = values$violinPlot, heatmap = values$heatmap,
      featurePlot_1 = values$featurePlots[["featurePlot_1"]], violinPlotGene_1 = values$violinPlotGenes[["violinPlotGene_1"]],
      featurePlot_2 = values$featurePlots[["featurePlot_2"]], violinPlotGene_2 = values$violinPlotGenes[["violinPlotGene_2"]],
      featurePlot_3 = values$featurePlots[["featurePlot_3"]], violinPlotGene_3 = values$violinPlotGenes[["violinPlotGene_3"]],
      featurePlot_4 = values$featurePlots[["featurePlot_4"]], violinPlotGene_4 = values$violinPlotGenes[["violinPlotGene_4"]],
      geneViolin = values$violin, annotation_umap = values$umap_annotation,
      sankey = values$sankey, plotmds = values$mds
    )
    values$saved_list[[new_index]] <- saved_list_tmp
  })

  # choose annotation colum from meta data
  observeEvent(input$annotation_column, {
    # print(input$annotation_column)
    values$annotations <- input$annotation_column
  })


  #feature plot and violin plot for gene
  observeEvent(input$genes_1, {
    values$featurePlots[["featurePlot_1"]] <- renderPlot(create_feature_plot(values$obj, input$genes_1, values))
    output$featurePlot_1 <- values$featurePlots[["featurePlot_1"]]
    values$violinPlotGenes[["violinPlotGene_1"]] <- renderPlot(create_violin_plot(values$obj, input$genes_1, values, ncol = NULL, pt.size = 0))
    output$violinPlotGene_1 <- values$violinPlotGenes[["violinPlotGene_1"]]
    values$selected_genes[["genes_1"]] = input$genes_1
  })
  observeEvent(input$genes_2, {
    values$featurePlots[["featurePlot_2"]] <- renderPlot(create_feature_plot(values$obj, input$genes_2, values))
    output$featurePlot_2 <- values$featurePlots[["featurePlot_2"]]
    values$violinPlotGenes[["violinPlotGene_2"]] <- renderPlot(create_violin_plot(values$obj, input$genes_2, values, ncol = NULL, pt.size = 0))
    output$violinPlotGene_2 <- values$violinPlotGenes[["violinPlotGene_2"]]
    values$selected_genes[["genes_2"]] = input$genes_2
  })
  observeEvent(input$genes_3, {
    values$featurePlots[["featurePlot_3"]] <- renderPlot(create_feature_plot(values$obj, input$genes_3, values))
    output$featurePlot_3 <- values$featurePlots[["featurePlot_3"]]
    values$violinPlotGenes[["violinPlotGene_3"]] <- renderPlot(create_violin_plot(values$obj, input$genes_3, values, ncol = NULL, pt.size = 0))
    output$violinPlotGene_3 <- values$violinPlotGenes[["violinPlotGene_3"]]
    values$selected_genes[["genes_3"]] = input$genes_3
  })
  observeEvent(input$genes_4, {
    values$featurePlots[["featurePlot_4"]] <- renderPlot(create_feature_plot(values$obj, input$genes_4, values))
    output$featurePlot_4 <- values$featurePlots[["featurePlot_4"]]
    values$violinPlotGenes[["violinPlotGene_4"]] <- renderPlot(create_violin_plot(values$obj, input$genes_4, values, ncol = NULL, pt.size = 0))
    output$violinPlotGene_4 <- values$violinPlotGenes[["violinPlotGene_4"]]
    values$selected_genes[["genes_4"]] = input$genes_4
  })

  # observeEvent(input$annotate, {
  #   show_modal_spinner(text = "updating annotations")
  #   clusters_list <- strsplit(input$clusters, "\\s+")
  #   if (input$annotation %in% values$annotations) {
  #     for (cluster in clusters_list) {
  #       values$annotations[[input$annotation]] <- append(values$annotations[[input$annotation]], cluster)
  #     }
  #   }
  #   else {
  #     values$annotations[[input$annotation]] = list()
  #     for (cluster in clusters_list) {
  #       values$annotations[[input$annotation]] <- append(values$annotations[[input$annotation]], cluster)
  #     }
  #   }
  #   # for (cluster in clusters_list) {
  #   #   values$annotation_show[[cluster]] <- input$annotation
  #   # }
  # })

  observeEvent(values$obj, {
    if (!is.null(values$obj)) {
      updateSelectizeInput(session, "genes_1", choices = rownames(values$obj), selected = values$selected_genes[["genes_1"]], server = TRUE)
      updateSelectizeInput(session, "genes_2", choices = rownames(values$obj), selected = values$selected_genes[["genes_2"]], server = TRUE)
      updateSelectizeInput(session, "genes_3", choices = rownames(values$obj), selected = values$selected_genes[["genes_3"]], server = TRUE)
      updateSelectizeInput(session, "genes_4", choices = rownames(values$obj), selected = values$selected_genes[["genes_4"]], server = TRUE)
      updateSelectizeInput(session, "annotation_column", choices = colnames(values$obj@meta.data))
      updateSelectInput(session, "clusters", choices = names(values$cluster_cell_counts))
    }
  })
}

# Run the application
shinyApp(ui, server)
