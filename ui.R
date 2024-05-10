source("functions/ClusteringFunc.R")
source("functions/InitialStepFunc.R")
source("tabs/ClusteringUI.R")
source("tabs/InitialStepUI.R")

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Single-Cell Clustering"),
  tabsetPanel(
    id = "main_tabs",
    InitialStepUI(),
    ClusteringUI()
  )
)
