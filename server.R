library(shiny)

if (!requireNamespace("shinyjs", quietly = TRUE)) install.packages("shinyjs")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("Biocmanager")
if (!requireNamespace(c("factoextra", "FactoMineR", "gProfileR", "PerformanceAnalytics"), quietly = TRUE)) install.packages(c("factoextra", "FactoMineR", "gProfileR", "PerformanceAnalytics"))
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("yaml", quietly = TRUE)) install.packages("yaml")
if (!requireNamespace("DT", quietly = TRUE)) install.packages("DT")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("UniProt.ws", quietly = TRUE)) BiocManager::install("UniProt.ws")
if (!requireNamespace("matrixTests", quietly = TRUE)) BiocManager::install("matrixTests")
if (!requireNamespace("tibble", quietly = TRUE)) BiocManager::install("tibble")


library(readr)
library(shinyjs)
library(yaml)
library(BiocManager)
library(artMS)
library(MSstats)
library(dplyr)
library(tidyr)
library(yaml)
library(DT)
library(UniProt.ws)
library(matrixTests)
library(tibble)


options(shiny.maxRequestSize = 1000*1024^2)

shinyServer(function(input, output, session) {
  
  observeEvent(input$generate_analysis, {
    req(input$evidencefile, input$keysfile, input$contrastfile)
    output$output_text <- renderText({"Working... Please Wait. (Might take a long time, be patient)"}) # Message for user patience 
    delay(1,{
    # Extract spectral counts
    spectral <- artmsSpectralCounts(
      evidence_file = input$evidencefile$datapath,
      keys_file = input$keysfile$datapath
    )
  
    # Read contrast file
    contrasts <- read_delim(input$contrastfile$datapath, delim = "-", col_names = c("condition1", "condition2"))
    
    # Initialize results list
    spectral <- spectral %>% filter(Proteins != "")
    results <- list()
    # Loop through each pair of conditions in the contrast file
    for (i in 1:nrow(contrasts)) {
      cond1 <- contrasts$condition1[i]
      cond2 <- contrasts$condition2[i]
     
      # Loop through each protein
      for (protein in unique(spectral$Proteins)) {
        data_cond1 <- spectral[spectral$Condition == cond1 & spectral$Proteins == protein, ]
        data_cond2 <- spectral[spectral$Condition == cond2 & spectral$Proteins == protein, ]
        
        # Check if there are enough observations for t-test
        if (nrow(data_cond1) > 1 & nrow(data_cond2) > 1) {
           t_test_results <- col_t_welch(data_cond1[,"spectral_counts"], data_cond2[,"spectral_counts"])
          
          
          # Add the result to the results list
          results[[length(results) + 1]] <- data.frame(
            Bait = cond1,
            Control= cond2,
            Prey = protein,
            SpecSum = paste(data_cond1$spectral_counts, collapse = "|"),
            ctrlCounts = paste(data_cond2$spectral_counts, collapse = "|"),
            p.value = t_test_results$pvalue
          )
        } 
      }
    }
    # Combine all results into a single data frame
    final_results <- do.call(rbind, results)
    
    # Save results to a text file
    write.table(final_results, file = "T-test.txt", sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Update UI with the table
    output$result_table <- renderDT({
      final_results
    }) 
    })
    
    # Enable download button
    output$download_results <- downloadHandler(
      filename = function() {
        "T-test.txt"
      },
      content = function(file) {
        write.table(final_results, file, sep = "\t", quote = FALSE, row.names = FALSE)
      }
    )
    
  })
  
})