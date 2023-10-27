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
      fileInput("data_file", "Choose h5ad File", accept = c(".h5ad")),
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
	 p("1. Upload gene expression matrix of myeloid cells",HTML("<br>"),"a. The matrix should be in h5ad format",HTML("<br>"),"b. The matrix should be normalized",HTML("<br>"),"c. You can view an example of generating a myeloid h5ad matrix by clicking",HTML("<a href='https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html'>here</a> and following the steps until you annotate the clusters. Then type the following:"),HTML("<br>"),"adata_myeloid = adata[adata.obs['leiden']=='myeloid',:]",HTML("<br>"),"adata.write_h5ad('Myeloid_Matrix.h5ad')"),
p("2. Once file is loaded, usage calculation will begin automatically"),
p("3. Once the 'Calculating Usages' process is completed, download the file by clicking 'Download Output Data'",HTML("<br>"),"    a. The output is automatically normalized by dividing each program usage score per cell by the sum of all usage scores for that cell and converting the values into percentages"),
p(),
h4("If you would like to run this program locally (may be needed for large matrices), you can download the program by", HTML("clicking <a href='https://github.com/BernsteinLab/Calculate_Myeloid_cNMF_Usage'>here</a> to visit the GitHub page of the tool.")),
p(),
h4("If you use this tool, please cite: Tyler E Miller, Chadi Abdul Kader El Farran, Charles P Couturier, et al., Programs, Origins, and Niches of Immunomodulatory Myeloid Cells in Gliomas. bioRxiv 2023.10.24.563466; doi:",  HTML("<a href='https://doi.org/10.1101/2023.10.24.563466'>https://doi.org/10.1101/2023.10.24.563466</a>"))
    )
  )
)

server <- function(input, output, session) {
  # Define Python environment
  reticulate::use_virtualenv("myenv")
  py <- import("your_python_script_module")  # Import your Python script module

  options(shiny.maxRequestSize = 3000 * 1024 ^ 2)

  observe({
    shinyFileChoose(input, "data_file", roots = c(base = ""))
  })

  # Function to process the data using Python
  observeEvent(input$data_file, {
  if (!is.null(input$data_file)) {
    data_file_path <- input$data_file$datapath

      # Call your Python function to process the H5AD data and get the Pandas DataFrame
      withProgress(message = 'Calculating Usages', value = 0, { processed_data <- py$process_data(data_file_path)

})
      # Create a downloadable link for the processed data
      output$download_output <- downloadHandler(
        filename = function() { paste("output_data", ".txt", sep = "\t") },
        content = function(file) {
          # Save the processed Pandas DataFrame as a CSV file
          write.table(processed_data, file, col.names=NA, sep = "\t", quote=FALSE)
        }
      )


    }
  })
}


shinyApp(ui, server)
