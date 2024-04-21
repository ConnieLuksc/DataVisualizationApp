source("global.R")

# Define UI
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Single-Cell Clustering"),
  tabsetPanel(
    id = "main_tabs",
    tabPanel(
      "Initial Step",
      fluidRow(
        column(
          2,
          helpText("Please upload RDS file without normalization and clustering"),
          fileInput("file", "Upload File", multiple = TRUE, accept = c(".rds")),
          fluidRow(
            column(6, numericInput("feature_upper", "nFeature_RNA", value = NA), style = "padding-right: 5px;"),
            column(6, numericInput("feature_lower", HTML("&nbsp;"), value = NA), style = "padding-left: 5px;")
          ),
          fluidRow(
            column(6, numericInput("count_upper", "nCount_RNA", value = NA), style = "padding-right: 5px;"),
            column(6, numericInput("count_lower", HTML("&nbsp;"), value = NA), style = "padding-left: 5px;")
          ),
          fluidRow(
            column(6, numericInput("percent_upper", "percent.mt", value = NA), style = "padding-right: 5px;"),
            column(6, numericInput("percent_lower", HTML("&nbsp;"), value = NA), style = "padding-left: 5px;")
          ),
          actionButton("filter_reset", "Reset"),
          actionButton("filter", "Filter"),
          selectInput("normalization_method", "Normalization Methods", c("LogNormalize", "CLR", "RC", "sctransform")),
          conditionalPanel(
            condition = "input.normalization_method == 'LogNormalize' || input.normalization_method == 'RC'",
            numericInput("parameter", NULL, value = NA)
          ),
          actionButton("normalize", "Normalize")
        ),
        column(
          10,
          fluidRow(
            column(5, plotOutput("filter_violinPlot")),
            column(7, plotOutput("feature_scatter"))
          ),
          fluidRow(
            column(5, plotOutput("elbowPlot")),
            column(7, plotOutput("feature_selection"))
          )
        )
      )
    ),
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
          # selectInput('clusters', 'Clusters', choices = NULL, multiple = TRUE, selectize = TRUE),
          # textInput("annotation", "Annotation"),
          # selectizeInput("gene", "Genes", choices = NULL),
          # actionButton("annotate", "Annotate"),
          actionButton("save", "Save", class = "save-btn"),
          actionButton("update", "Update", class = "update-btn")
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
  )
)

