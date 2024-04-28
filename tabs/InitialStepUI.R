InitialStepUI <- function() {
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
          numericInput("parameter", "Normalization Parameter", value = NA)
        ),
        numericInput("num_pcs", "PC value:", value = 30, min = 1, max = 100, step = 1),
        numericInput("num_features", "Number of Variable Features:", value = 2000, min = 100, max = 5000, step = 100),
        actionButton("normalize", "Normalize"),
        downloadButton("download_initial", "Download")
      ),
      column(
        10,
        fluidRow(
          column(5, plotOutput("filter_violinPlot")),
          column(7, plotOutput("feature_scatter"))
        ),
        fluidRow(
          column(7, plotOutput("feature_selection")),
          column(5, plotOutput("elbowPlot"))
        )
      )
    )
  )
}
