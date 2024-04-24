library(shiny)
library(shinyFiles)
library(reticulate)

# Create a Python virtual environment and install required packages
reticulate::virtualenv_create("myenv")
reticulate::virtualenv_install("myenv", packages=c("scipy", "scikit-learn==1.0.2", "numpy", "scanpy", "pandas"))

ui <- fluidPage(
  titlePanel("Calculate usages of the consensus cNMF programs in glioma-associated myeloid cells, as defined in Programs, Origins, and Niches of Immunomodulatory Myeloid Cells in Gliomas"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("mode_selector", "Choose Mode", choices = c("Annotation Mode", "Myeloid Program Calculation Mode")),
      selectInput("input_format", "Choose Input Format", choices = c("h5ad", "csv")),
      fileInput("data_file", "Choose File", accept = c(".h5ad", ".csv")),
      br(),
      downloadButton("download_output", "Download Output Data"),
      tags$style(HTML("
        p {
          margin-bottom: 40px
        }
      ")),
    ),
    mainPanel(
      h1("Welcome to the cNMF program usage calculator for glioma-associated myeloid cells"),
      h3("Instructions:"),
      p(),
      p("1. Choose the mode. Annotation Mode calculates enrichment of cell types NMF programs to help you annotate the Cells. Myeloid Program Calculation Mode calculates the usages of the consensus cNMF programs in glioma-associated myeloid cells. Upload gene expression matrix of myeloid cells for Myeloid Program Calculation Mode",HTML("<br>"),"a. Choose input format after choosing the mode. The matrix should be in h5ad or, csv. For CSV, genes in rows and cells in columns",HTML("<br>"),"b. The matrix should be normalized",HTML("<br>"),"c. You can view an example of generating a myeloid h5ad matrix by clicking",HTML("<a href='https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html'>here</a> and following the steps until you annotate the clusters. Then type the following:"),HTML("<br>"),"adata_myeloid = adata[adata.obs['leiden']=='myeloid',:]",HTML("<br>"),"adata.write_h5ad('Myeloid_Matrix.h5ad')"),
      p("2. Once the file is loaded, usage calculation will begin automatically"),
      p("3. Once the 'Calculating Usages' process is completed, download the file by clicking 'Download Output Data'",HTML("<br>"),"    a. The output is automatically normalized by dividing each program usage score per cell by the sum of all usage scores for that cell and converting the values into percentages"),
      p(),
      h4("If you would like to run this program locally (may be needed for large matrices), you can download the program by", HTML("clicking <a href='https://github.com/BernsteinLab/Calculate_Myeloid_cNMF_Usage'>here</a> to visit the GitHub page of the tool.")),
      p(),
      h4("This shiny app is developed and maintained by Chadi A. El Farran, M.Sc., Ph.D. (ChadiA_ElFarran@dfci.harvard.edu), with significant contributions from Charles P. Couturier, MD, Ph.D., and Tyler E. Miller, MD, Ph.D. If you use this tool, please cite: Tyler E Miller, Chadi Abdul Kader El Farran, Charles P Couturier, et al., Programs, Origins, and Niches of Immunomodulatory Myeloid Cells in Gliomas. bioRxiv 2023.10.24.563466; doi:",  HTML("<a href='https://doi.org/10.1101/2023.10.24.563466'>https://doi.org/10.1101/2023.10.24.563466</a>"))
    )
  )
)

server <- function(input, output, session) {
  # Define Python environment
  reticulate::use_virtualenv("myenv")
  py <- import("your_python_script_module")  # Import your Python script module

  options(shiny.maxRequestSize = 30000 * 1024 ^ 2)

  observe({
    shinyFileChoose(input, "data_file", roots = c(base = ""))
  })

  # Function to process the data using Python
  observeEvent(input$data_file, {
    if (!is.null(input$data_file)) {
      data_file_path <- input$data_file$datapath
      input_format <- input$input_format
      mode <- input$mode_selector

      # Call the appropriate Python script based on the selected mode
      withProgress(message = 'Calculating Usages', value = 0, {
        if (mode == "Annotation Mode") {
          if (input_format == "h5ad") {
            processed_data <- py$process_annotation_h5ad_data(data_file_path)
          } else if (input_format == "csv") {
            processed_data <- py$process_annotation_csv_data(data_file_path)
          } 
        } else if (mode == "Myeloid Program Calculation Mode") {
          if (input_format == "h5ad") {
            processed_data <- py$process_myeloid_program_h5ad_data(data_file_path)
          } else if (input_format == "csv") {
            processed_data <- py$process_myeloid_program_csv_data(data_file_path)
          } 
        }
      })

      # Create a downloadable link for the processed data
      output$download_output <- downloadHandler(
        filename = function() { paste("output_data", ".txt", sep = "\t") },
        content = function(file) {
          # Save the processed Pandas DataFrame as a CSV file
          write.table(processed_data, file, col.names=NA, sep = "\t", quote=FALSE, row.names=TRUE)
        }
      )
    }
  })
}

shinyApp(ui, server)
