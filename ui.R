library(shiny)
library(shinyjs)
library(DT)

shinyUI(fluidPage(
  useShinyjs(),
  titlePanel("T-test"),
  sidebarLayout(
    sidebarPanel(
      fileInput("evidencefile", "Upload Evidence File"),
      fileInput("keysfile", "Upload Keys File"),
      fileInput("contrastfile", "Upload Contrast File"),
      actionButton("generate_analysis", "Generate Analysis"),
      downloadButton("download_results", "Download Results")
    ),
    mainPanel(
      DTOutput("result_table")
    )
  )
))