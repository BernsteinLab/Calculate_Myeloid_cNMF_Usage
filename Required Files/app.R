library(shiny)
library(shinyFiles)
library(reticulate)

# Create a Python virtual environment and install required packages
reticulate::virtualenv_create("myenv")
reticulate::virtualenv_install("myenv", packages=c("scipy", "scikit-learn==1.0.2", "numpy", "scanpy", "pandas"))

ui <- fluidPage(
  titlePanel("Calculate usages of the consensus myeloid cNMF programs"),
  sidebarLayout(
    sidebarPanel(
      fileInput("data_file", "Choose h5ad File", accept = c(".h5ad")),
      br(),
      downloadButton("download_output", "Download Output Data")
    ),
    mainPanel(
	h1("Welcome to Myeloid-Glioma cNMF usage calculator"),
	h2("Instructions:"),
	p("1- The matrix should be in h5ad format"), 
p("2- The matrix should be normalized (It is recommended to perform sc.pp.normalize_per_cell and sc.pp.scale on a scanpy object)"),
p("3- Once the Calculating usages process is completed, you can proceed with downloading the file"),
p("4- The output is raw and should be normalized to percentages per cell to compare between cells"),
        h4("If you used this tool, please cite:")
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
