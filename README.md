# Calculate_Myeloid_cNMF_Usage
This shiny app receives a normalized h5ad matrix of myeloid cells and calculates the usage of the consensus myeloid cNMF programs in gliomas by the cells in the uploaded h5ad matrix.


The online version of this tool is available at: 
[(https://consensus-myeloid-program-calculator.shinyapps.io/shinyapp/)]

However, it is recommended that you download this app to avoid server delays and if you have large matrices.

The use of this tool is simple: Make sure to install all the required dependencies and follow the instructions below.

## Dependencies:
R > 4.0 (The latest version of RStudio is highly recommended)

Python > 3.7

The following R Libraries:

shiny

shinyFiles

reticulate

## App Setup:

1- Install the required R libraries (Needs to be done one time in each system):

>install.packages("shiny")

>install.packages("shinyFiles")

>install.packages("reticulate")

>reticulate::install_python(version = '3.9')


2- Download the four required files (**"app.R"**, **"Myeloid_NMF_Average_Gene_Spectra.txt"**, **"cnmf_run.spectra.k_18.dt_0_015.consensus.df.npz"**, and **"your_python_script_module.py"**) (Located in the "Required Files" folder) and **place them together in one folder** (This folder should be your working directory in R when you want to use the app). Ensure you download the files correctly through Git Hub (DO NOT right-click on the file to download it).   


3- Start your R session and set the working directory to the folder mentioned above in step 2 using setwd().


4- Load the required libraries as follows:

>library(shiny)

>library(shinyFiles)


5- Then type the following in R:

>shiny::runApp()


6- An interactive toolbox will be generated if everything is set up well.

7- Choose Annotation mode or Myeloid Calculation mode.

8— Upload the h5ad, csv or mtx matrix by clicking on browse and waiting for the "Calculating usage" bar to complete. (It should take a few minutes, depending on the size of the matrix.). If you have chosen mtx format, you will need to upload the features and barcodes files which are additional outputs of CellRanger or STARSolo.


9- Click "Download Output Data" to download the calculated usage.

## Input:

1- Choose the mode. "Annotation Mode" calculates the enrichment of cell types NMF programs to help you annotate the Cells. Myeloid Program Calculation Mode calculates the usages of the consensus cNMF programs in glioma-associated myeloid cells. Upload gene expression matrix of myeloid cells for Myeloid Program Calculation Mode.
a. Choose the input format after choosing the mode. The matrix should be in h5ad or csv. For CSV, genes should be in rows and cells in columns.
b. The matrix can be normalized or raw.
c. mtx and associated files can be gzipped or uncompressed.

2- Scanpy generated h5ad is highly recommended.

3- The values in adata.X can be normalized (It is recommended to use sc.pp.normalize_per_cell) or raw.

4- Genes should be stored as official gene symbols (HGNC official gene symbols).

Details are provided in the interactive toolbox.


## Output:

1-The output is a data frame in which the rows denote cells existing in the input h5ad matrix, whereas the columns represent the indicated consensus myeloid cNMF programs in human gliomas.

**2- IT IS IMPORTANT TO NOTE THAT THE USAGE VALUES IN THE OUTPUT ARE NORMALIZED AS PERCENTAGES PER CELL TO ENABLE COMPARISON AMONG CELLS**

Details are provided in the interactive toolbox.


## Demo Data:

Instructions for using the demo data are included in the "Demo Data" folder.

## Contact:

This shiny app is developed and maintained by Chadi A. El Farran, M.Sc., Ph.D. (ChadiA_ElFarran@dfci.harvard.edu), with significant contributions from Charles P. Couturier, MD, Ph.D., and Tyler E. Miller, MD, Ph.D.

If you use this tool, please cite:
