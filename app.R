
source('global.R')

ui <- dashboardPage(
    dashboardHeader(title = "scRNAseq Analysis"),
    dashboardSidebar(
        tags$head(
        tags$style(HTML(".skin-blue .main-header .sidebar-toggle {display: none;}"))
        ),
        sidebarMenu(id='tab',
            useShinyjs(),
            menuItem("Home Page", tabName = "home", icon = icon("list")),
            menuItem("scRNAseq Analyzer", tabName = "input", icon = icon("edit")),
            conditionalPanel(condition = "input.tab == 'input'",
                div(
                    fileInput("file", "Upload File", multiple=TRUE, accept=c('.rds')),
                    actionButton("reset", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                    actionButton("run", "Run", icon = icon("play"), style = "color: #fff; background-color: #28a745; width: 87.25%"),
                    numericInput("pc", "PC", value = NA),
                    numericInput("resolution", "Resolution", value = NA, step = 0.1),
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
                tags$h1(HTML("<u>Welcome to The scRNAseq Suerat analysis RShiny app</u>")),
                )
                )
            )         
)

server <- function(input, output, session) {
    options(shiny.maxRequestSize=300*1024^2)

    values <- reactiveValues(obj = NULL)

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

    observeEvent(input$run, {
        shinyjs::disable("run")

        # Clear tabs before 'Run' is ran another time
        removeTab("main_tabs", "UMAP")
        removeTab("main_tabs", "Gene Expression")
        removeTab("main_tabs", "Violin Plot")

        show_modal_spinner(text = "Preparing plots...")

        obj <- load_seurat_obj(input$file$datapath)
        values$obj <- obj
        if (is.vector(obj)){
            showModal(modalDialog(
                title = "Error with file",
                HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
                    paste(unlist(obj), collapse = "<br><br>"))
            ))
            shinyjs::enable("run")
            

        } else {
            
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
                    create_feature_plot(obj, input$gene)
                }
            })

            output$violinPlotGene <- renderPlot({
                if (!is.null(input$geneViolin)) {
                    create_violin_plot(obj, input$geneViolin)
                }
            })

            output$downloadFeaturePlot <- downloadHandler(
                filename = function(){
                    paste0(input$gene, '_feature_plot', '.png')
                },
                content = function(file){
                    plot <- create_feature_plot(obj, input$gene)
                    ggsave(filename=file, width = 10, height = 5, type = "cairo")
                }
            )
            output$download_umap <- downloadHandler(
                filename = function(){
                    paste0(input$metadata_col, '_UMAP', '.png')
                },
                content = function(file){
                    plot <- create_metadata_UMAP(obj, input$metadata_col)
                    ggsave(filename=file, width = 10, height = 5, type = "cairo")
                }
            )

            insertTab(
                inputId = "main_tabs",
                tabPanel(
                    "UMAP",
                    tabPanel("UMAP Plot",
                        fluidRow(
                            column(
                                width = 8,
                                plotOutput(outputId = 'umap'),
                                downloadButton("download_umap", "Download UMAP")
                            ),
                            # column(
                            #     width = 4,
                            #     selectizeInput("metadata_col", 
                            #         "Metadata Column", 
                            #         colnames(obj@meta.data)
                            #     )
                            # )
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
                        downloadButton("downloadFeaturePlot", "Download Feature Plot")
                    ),
                    column(
                        width = 4,
                        selectizeInput("gene", 
                            "Genes", 
                            rownames(obj)
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
                        downloadButton("downloadViolinPlot", "Download Violin Plot")
                    ),
                    column(
                        width = 4,
                        selectizeInput("geneViolin", 
                            "Genes", 
                            rownames(obj)
                        )
                    )
                    ),
                    style = "height: 90%; width: 95%; padding-top: 5%;"
                )
            )

            remove_modal_spinner()
            shinyjs::enable("run")
            
        }
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

}

shinyApp(ui, server)

