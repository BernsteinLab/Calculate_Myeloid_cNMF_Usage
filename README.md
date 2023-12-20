# Calculate_Myeloid_cNMF_Usage
This shiny app receives a normalized h5ad matrix of myeloid cells and calculates the usage of the consensus myeloid cNMF programs in gliomas by the cells in the uploaded h5ad matrix.


The online version of this tool is available at: 
[(https://consensus-myeloid-program-calculator.shinyapps.io/shinyapp/)]

However, it is recommended to download this app to avoid server delays and if you have large matrices.

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


2- Download the three required files (**"app.R"**, **"Myeloid_NMF_Average_Gene_Spectra.txt"**, and **"your_python_script_module.py"**) (Located in the "Required Files" folder) and **place them together in one folder** (This folder should be your working directory in R when you want to use the app).  Make sure you download the files properly through github (DO NOT right click on the file to download it).   


3- Start your R session and set the working directory to the folder mentioned above in step 2 using setwd().


4- Load the required libraries as follows:

>library(shiny)

>library(shinyFiles)


5- Then type the following in R:

>shiny::runApp()


6- An interactive toolbox will be generated if everything is set up well.


7- Upload the normalized h5ad matrix by clicking on browse and wait for the "Calculating usage" bar to complete. (It should take a few minutes, depending on the size of the matrix).


8- Then click on "Download Output Data" to download the calculated usages.

## Input:

1- Make sure that the cells in the input are myeloid in nature.

2- Scanpy generated h5ad is highly recommended.

3- The values in adata.X should be normalized (It is recommended to use sc.pp.normalize_per_cell).

4- Genes should be stored as official gene symbols (HGNC official gene symbols).

Details are provided in the interactive toolbox


## Output:

1-The output is a data frame in which the rows denote cells existing in the input h5ad matrix, whereas the columns represent the indicated consensus myeloid cNMF programs in human gliomas.

**2- IT IS IMPORTANT TO NOTE THAT THE USAGE VALUES IN THE OUTPUT ARE NORMALIZED AS PERCENTAGES PER CELL TO ENABLE COMPARISON AMONG CELLS**

Details are provided in the interactive toolbox


## Demo Data:

Instructions for using the demo data are included inside the "Demo Data" folder.

## Contact:

This shiny app is developed and maintained by Chadi A. El Farran, Ph.D. (ChadiA_ElFarran@dfci.harvard.edu), with significant contributions from Charles P. Couturier, MD, Ph.D., and Tyler E. Miller, MD, Ph.D.

If you use this tool, please cite:




