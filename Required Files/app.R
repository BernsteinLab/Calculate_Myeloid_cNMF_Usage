library(shiny)
library(reticulate)

# Use the virtual environment
use_virtualenv("myenv", required = TRUE)
py <- import("your_python_script_module")  # Ensure this matches your script name

ui <- fluidPage(
  titlePanel("cNMF Program Usage Calculator"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("mode_selector", "Choose Mode", choices = c("Annotation Mode", "Myeloid Program Calculation Mode")),
      selectInput("input_format", "Choose Input Format", choices = c("h5ad", "csv", "mtx")),
      uiOutput("file_input_ui"),
      br(),
      downloadButton("download_output", "Download Output Data"),
      tags$style(HTML("
        p {
          margin-bottom: 40px
        }
      "))
    ),
    mainPanel(
      h1("Welcome to the cNMF program usage calculator for human glioma expression data"),
      h3("Instructions:"),
      p(),
      p("1. Choose the mode. Annotation Mode calculates enrichment of cell types NMF programs to help you annotate the Cells. Myeloid Program Calculation Mode calculates the usages of the consensus cNMF programs in glioma-associated myeloid cells. Upload gene expression matrix of myeloid cells for Myeloid Program Calculation Mode",
        HTML("<br>"), "a. Choose input format after choosing the mode. The matrix should be in h5ad, csv, or mtx. For CSV, genes should be in rows and cells in columns. For mtx, you should upload the associated barcodes and features files which are outputs of CellRanger or STARSolo.",
        HTML("<br>"), "b. The matrix can be normalized or raw. mtx and associated files can be gz compressed or uncompressed",
        HTML("<br>"), "c. You can view an example of generating a myeloid h5ad matrix by clicking",
        HTML("<a href='https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html'>here</a> and following the steps until you annotate the clusters. Then type the following:"),
        HTML("<br>"), "adata_myeloid = adata[adata.obs['leiden']=='myeloid',:]",
        HTML("<br>"), "adata.write_h5ad('Myeloid_Matrix.h5ad')"),
      p("2. Once the file is loaded, usage calculation will begin automatically"),
      p("3. Once the 'Calculating Usages' process is completed, download the file by clicking 'Download Output Data'",
        HTML("<br>"), "a. The output is automatically normalized by dividing each program usage score per cell by the sum of all usage scores for that cell and converting the values into percentages"),
      p(),
      h4("If you would like to run this program locally (may be needed for large matrices), you can download the program by",
        HTML("clicking <a href='https://github.com/BernsteinLab/Calculate_Myeloid_cNMF_Usage'>here</a> to visit the GitHub page of the tool.")),
      p(),
      h4("This shiny app is developed and maintained by Chadi A. El Farran, M.Sc., Ph.D. (ChadiA_ElFarran@dfci.harvard.edu), with significant contributions from Charles P. Couturier, MD, Ph.D., and Tyler E. Miller, MD, Ph.D.",
        HTML("<br>"), "If you use this tool, please cite: Tyler E Miller, Chadi Abdul Kader El Farran, Charles P Couturier, et al., Programs, Origins, and Niches of Immunomodulatory Myeloid Cells in Gliomas. bioRxiv 2023.10.24.563466; doi:",
        HTML("<a href='https://doi.org/10.1101/2023.10.24.563466'>https://doi.org/10.1101/2023.10.24.563466</a>"))
    )
  )
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30000 * 1024 ^ 2)

  output$file_input_ui <- renderUI({
    if (input$input_format == "mtx") {
      tagList(
        fileInput("mtx_file", "Choose mtx File", accept = c(".mtx", ".mtx.gz")),
        fileInput("barcodes_file", "Choose Barcodes File", accept = c(".tsv", ".tsv.gz", ".txt", ".txt.gz")),
        fileInput("features_file", "Choose Features File", accept = c(".tsv", ".tsv.gz", ".txt", ".txt.gz"))
      )
    } else {
      fileInput("data_file", "Choose File", accept = c(".h5ad", ".csv"))
    }
  })

  observe({
    if (input$input_format == "mtx" && !is.null(input$mtx_file) && !is.null(input$barcodes_file) && !is.null(input$features_file)) {
      mtx_file_path <- input$mtx_file$datapath
      barcodes_file_path <- input$barcodes_file$datapath
      features_file_path <- input$features_file$datapath
      
      print(paste("MTX file path:", mtx_file_path))
      print(paste("Barcodes file path:", barcodes_file_path))
      print(paste("Features file path:", features_file_path))

      input_format <- input$input_format
      mode <- input$mode_selector

      print(paste("Input format:", input_format))
      print(paste("Mode:", mode))

      if (mode == "Annotation Mode") {
        print("Calling process_annotation_mtx_data")
        processed_data <- tryCatch({
          py$process_annotation_mtx_data(mtx_file_path, barcodes_file_path, features_file_path)
        }, error = function(e) {
          print(paste("Error in process_annotation_mtx_data:", e))
          return(NULL)
        })
      } else if (mode == "Myeloid Program Calculation Mode") {
        print("Calling process_myeloid_program_mtx_data")
        processed_data <- tryCatch({
          py$process_myeloid_program_mtx_data(mtx_file_path, barcodes_file_path, features_file_path)
        }, error = function(e) {
          print(paste("Error in process_myeloid_program_mtx_data:", e))
          return(NULL)
        })
      }

      print("Data processing completed for MTX format.")
      print(head(processed_data))

      output$data_output <- renderTable({
        head(processed_data)
      })

      output$download_output <- downloadHandler(
        filename = function() {
          paste("output_data", ".txt", sep = "")
        },
        content = function(file) {
          write.table(processed_data, file, row.names = TRUE, sep="\t", quote=FALSE, col.names=NA)
        }
      )
    } else if (!is.null(input$data_file)) {
      data_file_path <- input$data_file$datapath
      print(paste("Data file path:", data_file_path))

      input_format <- input$input_format
      mode <- input$mode_selector

      print(paste("Input format:", input_format))
      print(paste("Mode:", mode))

      if (mode == "Annotation Mode") {
        if (input_format == "h5ad") {
          print("Calling process_annotation_h5ad_data")
          processed_data <- tryCatch({
            py$process_annotation_h5ad_data(data_file_path)
          }, error = function(e) {
            print(paste("Error in process_annotation_h5ad_data:", e))
            return(NULL)
          })
        } else if (input_format == "csv") {
          print("Calling process_annotation_csv_data")
          processed_data <- tryCatch({
            py$process_annotation_csv_data(data_file_path)
          }, error = function(e) {
            print(paste("Error in process_annotation_csv_data:", e))
            return(NULL)
          })
        }
      } else if (mode == "Myeloid Program Calculation Mode") {
        if (input_format == "h5ad") {
          print("Calling process_myeloid_program_h5ad_data")
          processed_data <- tryCatch({
            py$process_myeloid_program_h5ad_data(data_file_path)
          }, error = function(e) {
            print(paste("Error in process_myeloid_program_h5ad_data:", e))
            return(NULL)
          })
        } else if (input_format == "csv") {
          print("Calling process_myeloid_program_csv_data")
          processed_data <- tryCatch({
            py$process_myeloid_program_csv_data(data_file_path)
          }, error = function(e) {
            print(paste("Error in process_myeloid_program_csv_data:", e))
            return(NULL)
          })
        }
      }

      print("Data processing completed for non-MTX format.")
      print(head(processed_data))

      output$data_output <- renderTable({
        head(processed_data)
      })

      output$download_output <- downloadHandler(
        filename = function() {
          paste("output_data", ".txt", sep = "")
        },
        content = function(file) {
          write.table(processed_data, file, row.names = TRUE, sep="\t", col.names=NA, quote=FALSE)
        }
      )
    } else {
      print("No files uploaded or files missing.")
      return(NULL)
    }
  })
}

shinyApp(ui, server)
