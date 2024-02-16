source("global.R")
library(shiny)
library(shinyjs)
library(DT)
library(purrr)
library(networkD3)

# Define UI
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Single-Cell Clustering"),
  tabsetPanel(
    id = "main_tabs",
    tabPanel(
      "Tab 1",
      fluidRow(
        column(
          2,
          fileInput("file", "Upload File", multiple = TRUE, accept = c(".rds")),
          actionButton("reset", "Reset", class = "reset-btn"),
          actionButton("run", "Run", class = "run-btn"),
          numericInput("pc", "PC", value = NA),
          numericInput("resolution", "Resolution", value = NA, step = 0.1),
          selectizeInput("gene", "Genes", choices = NULL),
          actionButton("save", "Save", class = "save-btn")
        ),
        column(
          8,
          fluidRow(
            column(5, plotOutput("violinPlot")),
            column(7, plotOutput("umap")),
          ),
          fluidRow(
            column(
              4,
              plotOutput(outputId = "featurePlot"),
            ),
            column(
              4,
              plotOutput(outputId = "violinPlotGene"),
            ),
            column(
              4,
              plotOutput(outputId = "mdsPlot")
            ),
            column(
              4,
              uiOutput(outputId = "sankeyPlot")
            ),
            column(
              4,
              plotOutput(outputId = "umap_annotation")
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
  values$selected_gene <- NULL


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
    values$selected_gene <- rownames(obj)[1]

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

      output$featurePlot <- renderPlot({
        if (!is.null(input$gene)) {
          create_feature_plot(values$obj, input$gene, values)
        }
      })

      output$violinPlot <- renderPlot({
        if (!is.na(input$pc) &&
          !is.na(input$resolution) &&
          !is.null(values$obj)) {
          values$obj[["percent.mt"]] <- PercentageFeatureSet(values$obj, pattern = "^MT-") # nolint
          violinPlot <- VlnPlot(values$obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # nolint
          values$violinPlot <- violinPlot
          violinPlot
        }
      })

      output$violinPlotGene <- renderPlot({
        if (!is.null(input$gene)) {
          create_violin_plot(values$obj, input$gene, values,
            ncol = NULL, pt.size = 0
          )
        }
      })

      output$cluster_cell_counts <- DT::renderDataTable(
        {
          cluster_ids <- Idents(values$obj)
          values$cluster_cell_counts <- table(cluster_ids)
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
        if (!is.na(input$pc) && !is.na(input$resolution)) {
          create_annotation_UMAP(
            obj, "seurat_clusters",
            input$pc, input$resolution, values
          )
        } else {
          ggplot() +
            theme_void() +
            geom_text(
              aes(
                x = 0.5, y = 0.5,
                label = str_wrap("Please input the parameters", width = 20)
              ),
              size = 12, color = "gray73", fontface = "bold"
            ) +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
        }
      })

      output$sankeyPlot <- renderUI({
        if (!is.null(values$obj)) {
          create_sankey_plot(values$obj, values, input$pc, input$resolution)
        }
      })


      output$text_output <- renderText({
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
              output[[paste0("featurePlot", values$count)]] <- renderPlot(current_saved_list[[key]]$feature) # nolint
              output[[paste0("geneViolin", values$count)]] <- renderPlot(current_saved_list[[key]]$geneViolin) # nolint
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
                          7,
                          plotOutput(paste0("violinPlot", values$count))
                        ),
                        column(
                          5,
                          plotOutput(paste0("umap", values$count))
                        ),
                      ),
                      fluidRow(
                        column(
                          6,
                          plotOutput(paste0("featurePlot", values$count)),
                        ),
                        column(
                          6,
                          plotOutput(paste0("geneViolin", values$count)),
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
      violin = values$violinPlot, feature = values$feature,
      geneViolin = values$violin, annotation_umap = values$umap_annotation,
      sankey = values$sankey
    )
    values$saved_list[[new_index]] <- saved_list_tmp
  })

  observeEvent(input$gene, {
    values$selected_gene <- input$gene
  })

  observeEvent(values$obj, {
    if (!is.null(values$obj)) {
      updateSelectizeInput(
        session, "gene",
        choices = rownames(values$obj),
        selected = values$selected_gene
      )
    }
  })
}

# Run the application
shinyApp(ui, server)
