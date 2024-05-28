ClusteringUI <- function() {
    tabPanel(
      "Clustering",
      fluidRow(
        column(
          2,
          #  fileInput("file", "Upload File", multiple = TRUE, accept = c('.rds')),
          actionButton("reset", "Reset", class = "reset-btn"),
          actionButton("run", "Run", class = "run-btn"),
          helpText("Please input resolution and PC before running:"),
          numericInput("pc", "PC", value = NA),
          numericInput("resolution", "Resolution", value = NA, step = 0.1),
          selectizeInput("genes_1", "Genes", choices = NULL, multiple = TRUE),
          selectizeInput("genes_2", NULL, choices = NULL, multiple = TRUE),
          selectizeInput("genes_3", NULL, choices = NULL, multiple = TRUE),
          selectizeInput("genes_4", NULL, choices = NULL, multiple = TRUE),
          numericInput("maximum_categories", label = "Maximum Categories", value = 100),
          selectizeInput("annotation_column", label = "annotation_column", choices = NULL),
          actionButton("save", "Save", class = "save-btn"),
          actionButton("update", "Update", class = "update-btn"),
          downloadButton("download_clustering", "Download"),
          downloadButton("download_r", "Download R object"),
          verbatimTextOutput("normalization"),
          verbatimTextOutput("normalization_parameter"),
          # verbatimTextOutput("PC_value"),
          verbatimTextOutput("Number_of_Variable_Features")

        ),
        column(
          8,
          fluidRow(
            column(5, plotOutput("violinPlot")),
            column(7, plotOutput("umap")),
          ),
          fluidRow(
            column(3, plotOutput(outputId = "featurePlot_1", height = "225px")),
            column(3, plotOutput(outputId = "featurePlot_2", height = "225px")),
            column(3, plotOutput(outputId = "featurePlot_3", height = "225px")),
            column(3, plotOutput(outputId = "featurePlot_4", height = "225px"))
          ),
          fluidRow(
            column(3, plotOutput(outputId = "violinPlotGene_1", height = "225px")),
            column(3, plotOutput(outputId = "violinPlotGene_2", height = "225px")),
            column(3, plotOutput(outputId = "violinPlotGene_3", height = "225px")),
            column(3, plotOutput(outputId = "violinPlotGene_4", height = "225px"))
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
          DT::dataTableOutput("cluster_cell_counts"),
          uiOutput("dynamic_elements")
        )
      ),
    )
}