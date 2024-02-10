source('global.R')
library(shiny)
library(shinyjs)
library(DT)
library(purrr)

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
               actionButton("save", "Save", class = "save-btn")
        ),
        column(8,
               fluidRow(
                 column(7, plotOutput("violinPlot")),
                 column(5, plotOutput("umap")),
               ),
               fluidRow(
                 column(
                   6,
                   plotOutput(outputId = 'featurePlot'),
                   selectizeInput("gene", "Genes", choices = NULL)
                 ),
                 column(
                   6,
                   plotOutput(outputId = 'violinPlotGene'),
                   selectizeInput("geneViolin", "Genes", choices = NULL)
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

  updateUI <- function(enable = TRUE) {
    if (enable) {
      shinyjs::enable("run")
      shinyjs::enable("pc")
      shinyjs::enable("resolution")
    } else {
      shinyjs::disable("run")
      shinyjs::disable("pc")
      shinyjs::disable("resolution")
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
          create_metadata_UMAP(obj, "seurat_clusters", input$pc, input$resolution, values)
        }else {
          ggplot() +
            theme_void() +
            geom_text(aes(x = 0.5, y = 0.5, label = str_wrap("Please input the parameters", width = 20)), size = 12, color = "gray73", fontface = "bold") +
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
          values$obj[["percent.mt"]] <- PercentageFeatureSet(values$obj, pattern = "^MT-")
          violinPlot <- VlnPlot(values$obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
          values$violinPlot <- violinPlot
          violinPlot
        }
      })

      output$violinPlotGene <- renderPlot({
        if (!is.null(input$geneViolin)) {
          create_violin_plot(values$obj, input$geneViolin, values)
        }
      })

      # count cell/cluster
      # output$cluster_cell_counts <- renderText({
      #     cluster_ids <- Idents(values$obj)
      #     cluster_cell_counts <- table(cluster_ids)
      #     values$result_text <- paste(names(cluster_cell_counts), ": ", as.numeric(cluster_cell_counts), collapse = ", ")
      #     return (values$result_text)
      # })

      output$cluster_cell_counts <- DT::renderDataTable({
        cluster_ids <- Idents(values$obj)
        values$cluster_cell_counts <- table(cluster_ids)
        data.frame(Cluster = names(values$cluster_cell_counts), Count = as.numeric(values$cluster_cell_counts))
      }, options = list(paging = FALSE, rownames = FALSE))


      output$text_output <- renderText({
        return("Current #Cells/Cluster")
      })

      shinyjs::enable("run")
    }
  })

  # Clear all sidebar inputs when 'Reset' button is clicked
  observeEvent(input$reset, {
    # Reset file input
    shinyjs::reset("file")

    # Clear plots and disable the 'Run' button
    output$umap <- renderPlot({ NULL })
    output$featurePlot <- renderPlot({ NULL })
    output$violinPlot <- renderPlot({ NULL })
    output$violinPlotGene <- renderPlot({ NULL })

    shinyjs::disable("run")
  })


  # monitor change of saved_list
  observeEvent(values$saved_list, {
    if (!is.null(values$saved_list)) {
      current_saved_list <- values$saved_list

      # update saved area part
      output$dynamic_elements <- renderUI({
        dynamic_elements <- lapply(seq_along(current_saved_list), function(key) {
          fluidRow(
            column(
              width = 12,
              actionButton(inputId = paste0("button_", key), label = paste("PC: ", current_saved_list[[key]]$pc, " Resolution: ", current_saved_list[[key]]$resolution)),
              # verbatimTextOutput(outputId = paste0("verbatim_output_", key))
              DT::dataTableOutput(outputId = paste0("verbatim_output_", key))

            )
          )
        })
        do.call(tagList, dynamic_elements)
      })

      lapply(seq_along(current_saved_list), function(key) {
        # output[[paste0("verbatim_output_", key)]] <- renderPrint({
        #     cat("Saved #Cell/Cluster", key, ": ", current_saved_list[[key]]$cluster)
        # })

        output[[paste0("verbatim_output_", key)]] <- DT::renderDataTable({
          data.frame(Cluster = names(current_saved_list[[key]]$cluster), Count = as.numeric(current_saved_list[[key]]$cluster))
        }, options = list(paging = FALSE, rownames = FALSE))

        observeEvent(input[[paste0("button_", key)]], {
          print("button")
          print(ignore_button_clicked)
          if (!ignore_button_clicked) {
            print("loaded")
            loaded_seurat <- load_seurat_obj(current_saved_list[[key]]$file)
            if (!is.vector(loaded_seurat)) {
              # values$obj <- loaded_seurat
              # updateNumericInput(session, "pc", value = current_saved_list[[key]]$pc)
              # updateNumericInput(session, "resolution", value = current_saved_list[[key]]$resolution)
              output[[paste0("umap", values$count)]] <- renderPlot(current_saved_list[[key]]$umap)
              output[[paste0("violinPlot", values$count)]] <- renderPlot(current_saved_list[[key]]$violin)
              output[[paste0("featurePlot", values$count)]] <- renderPlot(current_saved_list[[key]]$feature)
              output[[paste0("geneViolin", values$count)]] <- renderPlot(current_saved_list[[key]]$geneViolin)
              output[[paste0("textoutput", values$count)]] <- renderText({ return("Current #Cells/Cluster") })
              output[[paste0("cluster_cell_counts", values$count)]] <- DT::renderDataTable({
                cluster_ids <- Idents(values$obj)
                counts <- table(cluster_ids)
                data.frame(Cluster = names(counts), Count = as.numeric(counts))
              }, options = list(paging = FALSE))
              print(is.null(current_saved_list[[key]]$umap))
              insertTab(
                inputId = "main_tabs",
                tabPanel(
                  paste("Tab", values$count),
                  fluidRow(
                    column(2,
                           numericInput("pc", "PC", value = current_saved_list[[key]]$pc),
                           numericInput("resolution", "Resolution", value = current_saved_list[[key]]$resolution, step = 0.1),
                    ),
                    column(8,
                           fluidRow(
                             column(7, plotOutput(paste0("violinPlot", values$count))),
                             column(5, plotOutput(paste0("umap", values$count))),
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
                      # verbatimTextOutput("cluster_cell_counts"),
                      DT::dataTableOutput(paste0("cluster_cell_counts", values$count)),
                      # uiOutput("dynamic_elements")
                    )

                  )
                ),
                select = TRUE,
                # closable = TRUE,
              )
              values$count <- values$count + 1
              # shinyjs::runjs('$("#run").click();')
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
    # saved_list_tmp <- list(file = input$file$datapath, pc = input$pc, resolution = input$resolution, cluster = values$result_text)
    # values&umap <- output$umap()
    saved_list_tmp <- list(file = input$file$datapath, pc = input$pc, resolution = input$resolution, cluster = values$cluster_cell_counts, umap = values$umap, violin = values$violinPlot, feature = values$feature, geneViolin = values$violin)
    values$saved_list[[new_index]] <- saved_list_tmp
    # print(values$saved_list)
  })
  observe({
    if (!is.null(values$obj)) {
      updateSelectizeInput(session, "gene", choices = rownames(values$obj))
      updateSelectizeInput(session, "geneViolin", choices = rownames(values$obj))
    }
  })
}

# Run the application
shinyApp(ui, server)
