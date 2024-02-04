source('global.R')
library(shinyjs)
library(DT)
library(purrr)

ui <- dashboardPage(
    dashboardHeader(title = "Single-Cell Clustering"),
    dashboardSidebar(
        tags$head(
        tags$style(
            HTML("
            .main-header{
                background-color: #ffe;
            }
            .sidebar-toggle {display: true;}
            .reset-btn {
            color: #fff;
            background-color: #ff3b30;
            border-color: #ff3b30;
            padding: 10px 15px;
            font-size: 14px;
            font-weight: bold;
            border-radius: 6px;
            transition: background-color 0.3s;
            width: 87.25%;
          }
          .reset-btn:hover {
            background-color: #ff453a;
            border-color: #ff453a;
          }
          .run-btn {
            color: #fff;
            background-color: #28a745;
            border-color: #28a745;
            padding: 10px 15px;
            font-size: 14px;
            font-weight: bold;
            border-radius: 6px;
            transition: background-color 0.3s;
            width: 87.25%;
          }
          .run-btn:hover {
            background-color: #218838;
            border-color: #218838;
          }
          .apple-download-btn {
          background-color: #4CAF50;
          border: none;
          color: white;
          text-align: center;
          text-decoration: none;
          display: inline-block;
          font-size: 16px;
          margin: 4px 2px;
          cursor: pointer;
          border-radius: 5px;
          transition: background-color 0.3s;
        }

        .apple-download-btn:hover {
          background-color: #45a049;
        }


            ")

        )
        ),
        sidebarMenu(id='tab',
            useShinyjs(),
            menuItem("Home", tabName = "home", icon = icon("list")),
            menuItem("Single-Cell Clustering Analyzer", tabName = "input", icon = icon("edit")),
            conditionalPanel(condition = "input.tab == 'input'",
                div(
                    # upload file
                    fileInput("file", "Upload File", multiple=TRUE, accept=c('.rds')),
                    actionButton("reset", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                    actionButton("run", "Run", icon = icon("play"), style = "color: #fff; background-color: #28a745; width: 87.25%"),
                    numericInput("pc", "PC", value = NA),
                    numericInput("resolution", "Resolution", value = NA, step = 0.1),
                    actionButton("save", "Save", style = "color: #fff; background-color: #dc3545; width: 87.25%")
                    )
                )
        )
    ),
    dashboardBody(
            tabItems(
                tabItem(tabName = "input", # tabItem refers to tab in sidebar (not main panel)
                tabsetPanel(id = 'main_tabs',
                    tabPanel("Instructions",
                            includeMarkdown("./markdown/instructions.md")
                            )
                        )
                    ),
                tabItem(tabName = "home",
                tags$h1(HTML("<u>Welcome to The Single-Cell Clustering Web App</u>"))
                )
                )
            ),
            skin = "purple"
)

server <- function(input, output, session) {
    options(shiny.maxRequestSize=300*1024^2)
    clicked_link <- reactiveVal(FALSE)
    values <- reactiveValues()
    saved_list <- list(
      # item1 = list(file = "/Users/pengruiyang/Desktop/Shiny/single_cell_shiny_app/pbmc_tutorial.rds", cluster = "0: 1, 1: 2"),
      # item2 = list(file = "/Users/pengruiyang/Desktop/Shiny/single_cell_shiny_app/pbmc_tutorial.rd", cluster = "0: 2, 1: 3")
    )
    # Disable Run by default
    shinyjs::disable("run")
    shinyjs::disable("pc")
    shinyjs::disable("resolution")

    observe({
    if(is.null(input$file) != TRUE) {
        shinyjs::enable("run")
        shinyjs::enable("pc")
        shinyjs::enable("resolution")
    } else {
        shinyjs::disable("run")
        shinyjs::disable("pc")
        shinyjs::disable("resolution")
    }
    })

    output$dynamic_elements <- renderUI({
      # current_saved_list <-saved_list()
      dynamic_elements <- imap(saved_list, function(item, key) {
        fluidRow(
          column(
            width = 12,
            actionButton(inputId = paste0("button_", key), label = paste("Button ", key)),
            verbatimTextOutput(outputId = paste0("verbatim_output_", key))
        )
        )
      })
      do.call(tagList, dynamic_elements)
    })
    # current_saved_list2 <- saved_list()
    lapply(names(saved_list), function(key) {
      output[[paste0("verbatim_output_", key)]] <- renderPrint({
        cat(" Saved #Cell/Cluster ", key, ": ", saved_list[[key]]$cluster)
      })
    })

    observeEvent(input$run, {
        shinyjs::disable("run")

        # Clear tabs before 'Run' is ran another time
        removeTab("main_tabs", "UMAP")
        removeTab("main_tabs", "Gene Expression")
        removeTab("main_tabs", "Violin Plot")

        show_modal_spinner(text = "Preparing plots...")

        obj <- load_seurat_obj(input$file$datapath)
        values$obj <- load_seurat_obj(input$file$datapath)
        values$run_triggered <- reactiveVal(FALSE)
        if (is.vector(values$obj)){
            showModal(modalDialog(
                title = "Error with file",
                HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
                    paste(unlist(values$obj), collapse = "<br><br>"))
            ))
            shinyjs::enable("run")


        } else {
            # umap
            output$umap <- renderPlot({
                if (!is.na(input$pc) && !is.na(input$resolution)) {
                    show_modal_spinner(text = "Preparing plots...")
                    create_metadata_UMAP(obj, "seurat_clusters", input$pc, input$resolution, values)
                }else{
                ggplot() +
                    theme_void() +
                    geom_text(aes(x = 0.5, y = 0.5, label = str_wrap("Please input the parameters", width = 20)), size = 12, color = "gray73", fontface = "bold") +
                    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
                }
            })

            # output$violinPlot <- renderPlot({
            #     if (!is.na(input$pc) && !is.na(input$resolution)) {
            #         values$obj[["percent.mt"]] <- PercentageFeatureSet(values$obj, pattern = "^MT-")
            #         VlnPlot(values$obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
            #     }
            # })

            output$featurePlot <- renderPlot({
                if (!is.null(input$gene)) {
                    create_feature_plot(values$obj, input$gene)
                }
            })

            output$violinPlotGene <- renderPlot({
                if (!is.null(input$geneViolin)) {
                    create_violin_plot(values$obj, input$geneViolin)
                }
            })

            output$downloadFeaturePlot <- downloadHandler(
                filename = function(){
                    paste0(input$gene, '_feature_plot', '.png')
                },
                content = function(file){
                    plot <- create_feature_plot(values$obj, input$gene)
                    ggsave(filename=file, width = 10, height = 5, type = "cairo")
                }
            )
            output$download_umap <- downloadHandler(
                filename = function(){
                    paste0(input$metadata_col, '_UMAP', '.png')
                },
                content = function(file){
                    plot <- create_metadata_UMAP(values$obj, input$metadata_col)
                    ggsave(filename=file, width = 10, height = 5, type = "cairo")
                }
            )
          # count cell/cluster
            output$cluster_cell_counts <- renderText({
              cluster_ids <- Idents(values$obj)
              cluster_cell_counts <- table(cluster_ids)
              values$result_text <- paste(names(cluster_cell_counts), ": ", as.numeric(cluster_cell_counts), collapse = ", ")
              return (values$result_text)
              })

          output$text_output <- renderText({
            return("Current #Cells/Cluster")
          })


            insertTab(
                inputId = "main_tabs",
                tabPanel(
                    "UMAP",
                    tabPanel("UMAP Plot",
                        fluidRow(
                            column(
                                width = 4,
                                plotOutput(outputId = 'umap'),
                                downloadButton("download_umap", "Download UMAP", class = "apple-download-btn")
                            ),
                            column(
                                width = 4,
                                selectizeInput("metadata_col",
                                    "Metadata Column",
                                    colnames(values$obj@meta.data),
                                )
                            ),
                            # saved area
                            column(
                                width = 4,
                                style = "overflow-y: scroll; max-height: 300px;",
                                textOutput(outputId = "text_output"),
                                verbatimTextOutput("cluster_cell_counts"),
                                uiOutput("dynamic_elements")

                            )
                        ),
                        style = "height: 90%; width: 95%; padding-top: 5%;"
                    ),
                    tabPanel("Violin Plot",
                        fluidRow(
                            column(
                                width = 12,
                                plotOutput(outputId = 'violinPlot')  # Create an output for the violin plot
                            )
                        ),
                        style = "height: 90%; width: 95%; padding-top: 5%;"
                    )

                ),
                select = TRUE
            )
            insertTab(
                inputId = "main_tabs",
                tabPanel(
                    "Gene Expression",
                    fluidRow(
                    column(
                        width = 8,
                        plotOutput(outputId = 'featurePlot'),
                        downloadButton("downloadFeaturePlot", "Download Feature Plot", class = "apple-download-btn")
                    ),
                    column(
                        width = 4,
                        selectizeInput("gene",
                            "Genes",
                            rownames(values$obj)
                        )
                    )
                    ),
                    style = "height: 90%; width: 95%; padding-top: 5%;"
                )
            )
            insertTab(
                inputId = "main_tabs",
                tabPanel(
                    "Violin Plot",
                    fluidRow(
                    column(
                        width = 8,
                        plotOutput(outputId = 'violinPlotGene'),
                        downloadButton("downloadViolinPlot", "Download Violin Plot", class = "apple-download-btn")
                    ),
                    column(
                        width = 4,
                        selectizeInput("geneViolin",
                            "Genes",
                            rownames(values$obj)
                        )
                    )
                    ),
                    style = "height: 90%; width: 95%; padding-top: 5%;"
                )
            )

            remove_modal_spinner()
            shinyjs::enable("run")

            clicked_link(FALSE)

        }
    })


  # current_saved_list1 <- saved_list()
  # build observeEvent for the buttons in the list
  lapply(names(saved_list), function(key) {
    observeEvent(input[[paste0("button_", key)]], {
      isolate({
        print(paste("button_", key))

      local_key <- key
      if (!clicked_link()) {
        shinyjs::disable("reset-btn")

        # Extract file path from saved_list
        file_path <- saved_list[[local_key]]$file

        # Load Seurat object based on the file_path
        loaded_seurat <- load_seurat_obj(file_path)

        # Check if loading is successful
        if (is.vector(loaded_seurat)) {
          showModal(modalDialog(
            title = "Error with file",
            HTML(paste("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
                       paste(unlist(loaded_seurat), collapse = "<br><br>"))
          )))
          shinyjs::enable("reset-btn")
        } else {
          # Update values$obj
          values$obj <- loaded_seurat
          clicked_link(TRUE)
          shinyjs::runjs('$("#run").click();')
        }
      }
      })
    })
})

    observeEvent(values$obj, {
            output$violinPlot <- renderPlot({
                if (!is.na(input$pc) && !is.na(input$resolution) && !is.null(values$obj)) {
                    values$obj[["percent.mt"]] <- PercentageFeatureSet(values$obj, pattern = "^MT-")
                    VlnPlot(values$obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
                }
            })

    })

    # Clear all sidebar inputs when 'Reset' button is clicked
    observeEvent(input$reset, {
        shinyjs::reset("file")
        removeTab("main_tabs", "UMAP")
        removeTab("main_tabs", "Gene Expression")
        removeTab("main_tabs", "Violin Plot")
        shinyjs::disable("run")
    })

    observeEvent(input$save, {
      print(saved_list)
      observe({
  cat("Debugging - Length of saved_list: ", length(saved_list), "\n")
})
        new_index <- paste0("item", length(saved_list) + 1)
        saved_list_tmp <- list(file = input$file$datapath, pc = input$pc, resolution = input$resolution, cluster = values$result_text)
        saved_list[[new_index]] <<- saved_list_tmp
        print(saved_list)


        output$dynamic_elements <- renderUI({
            dynamic_elements <- imap(saved_list, function(item, key) {
            fluidRow(
                column(
                width = 12,
                actionButton(inputId = paste0("button_", key), label = paste("PC: ", saved_list[[key]]$pc, " Resolution: ", saved_list[[key]]$resolution)),
                verbatimTextOutput(outputId = paste0("verbatim_output_", key))
            )
            )
            })
            do.call(tagList, dynamic_elements)
        })
          lapply(names(saved_list), function(key) {
          output[[paste0("verbatim_output_", key)]] <- renderPrint({
            cat(" Saved #Cell/Cluster ", key, ": ", saved_list[[key]]$cluster)
          })
        })
          lapply(names(saved_list), function(key) {
            observeEvent(input[[paste0("button_", key)]], {
              isolate({
                print(paste("button_", key))

              local_key <- key
              if (!clicked_link()) {
                shinyjs::disable("reset-btn")

                # Extract file path from saved_list
                file_path <- saved_list[[local_key]]$file

                # Load Seurat object based on the file_path
                loaded_seurat <- load_seurat_obj(file_path)

                # Check if loading is successful
                if (is.vector(loaded_seurat)) {
                  showModal(modalDialog(
                    title = "Error with file",
                    HTML(paste("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
                               paste(unlist(loaded_seurat), collapse = "<br><br>"))
                  )))
                  shinyjs::enable("reset-btn")
                } else {
                  # Update values$obj
                  values$obj <- loaded_seurat
                  clicked_link(TRUE)
                  shinyjs::runjs('$("#run").click();')
                }
              }
              })
            })
    })
    })


}

shinyApp(ui, server)

