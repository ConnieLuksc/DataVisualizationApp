
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
  )
}
