library(shiny)
library(reticulate)


# Check if virtual environment exists, if not create it and install required packages
if (!reticulate::virtualenv_exists("myenv")) {
  reticulate::virtualenv_create("myenv")
  reticulate::virtualenv_install("myenv", packages = c("cnmf", "scanpy", "pandas", "numba"))

}

# Use the virtual environment
reticulate::use_virtualenv("myenv", required = TRUE)

# Ensure the Python script module is imported correctly
py <- NULL
try({
  py <- import("your_python_script_module", convert = TRUE)
}, silent = FALSE)

ui <- fluidPage(
  titlePanel("cNMF Program Usage Calculator"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("mode_selector", "Choose Mode", choices = c("Annotation Mode", "Myeloid Program Calculation Mode")),
      selectInput("input_format", "Choose Input Format", choices = c("h5ad", "csv", "mtx")),
      uiOutput("file_input_ui"),
      br(),
      downloadButton("download_output", "Download Output Data"),
      downloadButton("download_log", "Download Log File"),
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
      p("1. Choose the mode. Annotation Mode calculates the enrichment of cell types NMF programs to help you annotate the Cells. Myeloid Program Calculation Mode calculates the usages of the consensus cNMF programs in glioma-associated myeloid cells. Upload gene expression matrix of myeloid cells for Myeloid Program Calculation Mode",
        HTML("<br>"), "a. Choose input format after choosing the mode. The matrix should be in h5ad, csv, or mtx. For CSV, genes should be in rows and cells in columns. For mtx, you should upload the associated barcodes and feature files, which are outputs of CellRanger or STARSolo. Feature file has to have three columns, with the second column including gene symbols and the third column having the words Gene Expression in all rows",
        HTML("<br>"), "b. The matrix should be normalized. mtx and associated files can be gz compressed or uncompressed. To make the process faster, upload gzipped mtx file",
        HTML("<br>"), "c. You can view an example of generating a myeloid h5ad matrix by clicking",
        HTML("<a href='https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html'>here</a> and following the steps until you annotate the clusters. Then type the following:"),
        HTML("<br>"), "adata_myeloid = adata[adata.obs['leiden']=='myeloid',:]",
        HTML("<br>"), "adata.write_h5ad('Myeloid_Matrix.h5ad')"),
      p("2. Once the file(s) is/are loaded, usage calculation will begin automatically."),
      p("3. Once the 'Calculating Usages' process is completed, the progress bar will disappear, and you can download the output file by clicking 'Download Output Data'. In case of an error, you can download the log file which can help in troubleshooting.",
        HTML("<br>"), "a. The output is automatically normalized by dividing each program usage score per cell by the sum of all usage scores for that cell and converting the values into percentages"),
      p(),
      h4("If you would like to run this program locally (may be needed for large matrices), you can download the program by",
        HTML("clicking <a href='https://github.com/BernsteinLab/Calculate_Myeloid_cNMF_Usage'>here</a> to visit the GitHub page of the tool.")),
      p(),
      h4("This shiny app is developed and maintained by Chadi A. El Farran, M.Sc., Ph.D. (Chadi.ElFarran@stjude.org), with significant contributions from Charles P. Couturier, MD, Ph.D., and Tyler E. Miller, MD, Ph.D.",
        HTML("<br>"), "If you use this tool, please cite: Miller, T.E., El Farran, C.A., Couturier, C.P. et al. Programs, origins and immunomodulatory functions of myeloid cells in glioma. Nature (2025); doi:",
        HTML("<a href='https://doi.org/10.1038/s41586-025-08633-8'>https://doi.org/10.1038/s41586-025-08633-8</a>"))
    )
  )
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30000 * 1024 ^ 2)

Sys.setenv(LD_LIBRARY_PATH = "/opt/homebrew/opt/zlib/lib:/usr/local/lib")
Sys.setenv(DYLD_LIBRARY_PATH = "/opt/homebrew/opt/zlib/lib:/usr/local/lib")
  
  log_file <- tempfile(fileext = ".txt")
  
  log_message <- function(message) {
    if (is.list(message)) {
      cat(capture.output(str(message)), file = log_file, append = TRUE, sep = "\n")
    } else {
      cat(message, "\n", file = log_file, append = TRUE)
    }
  }
  
  output$file_input_ui <- renderUI({
    if (input$input_format == "mtx") {
      tagList(
        fileInput("mtx_file", "Choose mtx File", accept = c(".mtx", ".mtx.gz")),
        fileInput("barcodes_file", "Choose Barcodes File", accept = c(".tsv", ".gz", ".txt")),
        fileInput("features_file", "Choose Features File", accept = c(".tsv", ".gz", ".txt"))
      )
    } else {
      fileInput("data_file", "Choose File", accept = c(".h5ad", ".csv"))
    }
  })

  observe({
    withProgress(message = 'Processing', value = 0, {
      incProgress(0.1, detail = "Preparing files")
      
      if (input$input_format == "mtx" && !is.null(input$mtx_file) && !is.null(input$barcodes_file) && !is.null(input$features_file)) {
        log_message("MTX format detected")
        mtx_file_path <- input$mtx_file$datapath
        barcodes_file_path <- input$barcodes_file$datapath
        features_file_path <- input$features_file$datapath
        
        log_message(paste("MTX file path:", mtx_file_path))
        log_message(paste("Barcodes file path:", barcodes_file_path))
        log_message(paste("Features file path:", features_file_path))

        input_format <- input$input_format
        mode <- input$mode_selector

        log_message(paste("Input format:", input_format))
        log_message(paste("Mode:", mode))

        incProgress(0.4, detail = "Running analysis")

        if (mode == "Annotation Mode") {
          log_message("Calling process_annotation_mtx_data")
          processed_data <- tryCatch({
            py$process_annotation_mtx_data(mtx_file_path, barcodes_file_path, features_file_path)
          }, error = function(e) {
            log_message(paste("Error in process_annotation_mtx_data:", e))
            return(NULL)
          })
        } else if (mode == "Myeloid Program Calculation Mode") {
          log_message("Calling process_myeloid_program_mtx_data")
          processed_data <- tryCatch({
            py$process_myeloid_program_mtx_data(mtx_file_path, barcodes_file_path, features_file_path)
          }, error = function(e) {
            log_message(paste("Error in process_myeloid_program_mtx_data:", e))
            return(NULL)
          })
        }

        incProgress(0.9, detail = "Finalizing")

        if (!is.null(processed_data)) {
          log_message("Data processing completed for MTX format.")
          log_message(capture.output(head(processed_data)))

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
        }
        
      } else if (!is.null(input$data_file)) {
        log_message("Non-MTX format detected")
        data_file_path <- input$data_file$datapath
        log_message(paste("Data file path:", data_file_path))

        input_format <- input$input_format
        mode <- input$mode_selector

        log_message(paste("Input format:", input_format))
        log_message(paste("Mode:", mode))

        incProgress(0.4, detail = "Running analysis")

        if (mode == "Annotation Mode") {
          if (input_format == "h5ad") {
            log_message("Calling process_annotation_h5ad_data")
            processed_data <- tryCatch({
              py$process_annotation_h5ad_data(data_file_path)
            }, error = function(e) {
              log_message(paste("Error in process_annotation_h5ad_data:", e))
              return(NULL)
            })
          } else if (input_format == "csv") {
            log_message("Calling process_annotation_csv_data")
            processed_data <- tryCatch({
              py$process_annotation_csv_data(data_file_path)
            }, error = function(e) {
              log_message(paste("Error in process_annotation_csv_data:", e))
              return(NULL)
            })
          }
        } else if (mode == "Myeloid Program Calculation Mode") {
          if (input_format == "h5ad") {
            log_message("Calling process_myeloid_program_h5ad_data")
            processed_data <- tryCatch({
              py$process_myeloid_program_h5ad_data(data_file_path)
            }, error = function(e) {
              log_message(paste("Error in process_myeloid_program_h5ad_data:", e))
              return(NULL)
            })
          } else if (input_format == "csv") {
            log_message("Calling process_myeloid_program_csv_data")
            processed_data <- tryCatch({
              py$process_myeloid_program_csv_data(data_file_path)
            }, error = function(e) {
              log_message(paste("Error in process_myeloid_program_csv_data:", e))
              return(NULL)
            })
          }
        }

        incProgress(0.9, detail = "Finalizing")

        if (!is.null(processed_data)) {
          log_message("Data processing completed for non-MTX format.")
          log_message(capture.output(head(processed_data)))

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
        }
      } else {
        log_message("No files uploaded or files missing.")
        return(NULL)
      }
      
      incProgress(1, detail = "Done")
    })
  })
  
  output$download_log <- downloadHandler(
    filename = function() {
      "log_file.txt"
    },
    content = function(file) {
      file.copy(log_file, file)
    }
  )
}

shinyApp(ui, server)
