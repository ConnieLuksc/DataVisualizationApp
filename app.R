#version 1; 20240129
source('global.R')

callFB <- "
shinyjs.FB = function() {
    var firebaseConfig = {
        apiKey: 'AIzaSyDAEi8Jjq7mYo5ViHfELPR83GVV0F',
        authDomain: 'scrnaseq-app.firebaseapp.com',
        databaseURL: 'https://scrnaseq-app-default-rtdb.firebaseio.com',
        projectId: 'scrnaseq-app',
        storageBucket: 'scrnaseq-app.appspot.com',
        messagingSenderId: '670365424524',
        appId: '1:670365424524:web:40bf001b36b573ba0814a9',
        measurementId: 'G-HPR83GVV0F'
    };
    firebase.initializeApp(firebaseConfig);
    var database = firebase.database();
    var id = database.ref('test').push().key;
    var myData = 'hello';
    database.ref('test/' + id + '/').update({myData}).then(function() {
        console.log('myData sent!');
    });
}";


ui <- dashboardPage(
    dashboardHeader(title = "scRNAseq Analysis"),
    dashboardSidebar(
        tags$head(
        tags$style(HTML(".skin-blue .main-header .sidebar-toggle {display: none;}"))
        ),
        sidebarMenu(id='tab',
            useShinyjs(),
            extendShinyjs(text = callFB, functions = c("FB")),
            tags$head(tags$script(src = "https://www.gstatic.com/firebasejs/8.6.1/firebase-app.js")),
            tags$head(tags$script(src = "https://www.gstatic.com/firebasejs/8.6.1/firebase-database.js")),
            menuItem("Home Page", tabName = "home", icon = icon("list")),
            menuItem("scRNAseq Analyzer", tabName = "input", icon = icon("edit")),
            conditionalPanel(condition = "input.tab == 'input'",
                div(
                    fileInput("file", "Upload File", multiple = TRUE, accept = c('.rds')),
                    actionButton("reset", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                    actionButton("run", "Run", icon = icon("play"), style = "color: #fff; background-color: #28a745; width: 87.25%"),
                    actionButton("save_to_firebase", "Save to Firebase", icon = icon("cloud-upload"), style = "margin-top: 5px;"),
                    actionButton("load_from_firebase", "Load from Firebase", icon = icon("cloud-download"), style = "margin-top: 5px;")
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

    values <- reactiveValues()

    # Disable Run by default
    shinyjs::disable("run")

    observe({
    if(is.null(input$file) != TRUE) {
        shinyjs::enable("run")
    } else {
        shinyjs::disable("run")
    }
    })

    observeEvent(input$run, {
        shinyjs::disable("run")

        # Clear tabs before 'Run' is ran another time
        removeTab("main_tabs", "UMAP")
        removeTab("main_tabs", "Gene Expression")

        show_modal_spinner(text = "Preparing plots...")

        obj <- load_seurat_obj(input$file$datapath)
        if (is.vector(obj)){
            showModal(modalDialog(
                title = "Error with file",
                HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
                    paste(unlist(obj), collapse = "<br><br>"))
            ))
            shinyjs::enable("run")

        } else {
            
            output$umap <- renderPlot({
                if (!is.null(input$metadata_col)) {
                    create_metadata_UMAP(obj, input$metadata_col)
                }
            })
            output$violinPlot <- renderPlot({
                if (!is.null(obj)) {
                    VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
                }
            })

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
                            column(
                                width = 4,
                                selectizeInput("metadata_col", 
                                    "Metadata Column", 
                                    colnames(obj@meta.data)
                                )
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

    # Clear all sidebar inputs when 'Reset' button is clicked
    observeEvent(input$reset, {
        shinyjs::reset("file")
        removeTab("main_tabs", "UMAP")
        removeTab("main_tabs", "Gene Expression")
        shinyjs::disable("run")
    })

    observeEvent(input$save_to_firebase, {
        shinyjs::js$FB()  # Correct way to call the JavaScript function
    })

        # if (exists("obj")) {
        #     json_data <- toJSON(obj)
        #     firebase_url <- "https://scrnaseq-app-default-rtdb.firebaseio.com/data.json" # replace with your URL

        #     response <- PUT(
        #       url = firebase_url,
        #       body = json_data,
        #       encode = "json"
        #     )

        #     if (http_status(response)$category == "success") {
        #         showNotification("Data saved to Firebase successfully!", type = "message")
        #     } else {
        #         showNotification("Failed to save data to Firebase", type = "error")
        #     }
        # } else {
        #     showNotification("No data available to save", type = "error")
        # }
    # })

    observeEvent(input$load_from_firebase, {
        firebase_url <- "https://scrnaseq-app-default-rtdb.firebaseio.com/data.json" # replace with your URL

        response <- GET(firebase_url)

        if (http_status(response)$category == "success") {
            json_data <- content(response, "text")
            obj <- fromJSON(json_data)
            showNotification("Data loaded from Firebase successfully!", type = "message")
        } else {
            showNotification("Failed to load data from Firebase", type = "error")
        }
    })

}

shinyApp(ui, server)

