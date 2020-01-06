# **The aging skin microenvironment dictates stem cell behavior**

This project provides the code developed in the study of Yejing Ge _et al._ **_"The aging skin microenvironment dictates stem cell behavior"_**. It focuses on comparing Young vs. Aged mouse skin with single-cell RNA-seq.

## **Requirements**

* R (tested in R version 3.5.2 (2018-12-20) -- "Eggshell Igloo")
* R librarys: Seurat (v2.3), SingleR (github version)

## **Data**

Raw counts data and processed Seurat object will be released after the final publication.

## **Reproduce results**

#### **1. Data preprocess**
Move Cell ranger output foder under `data`.
[1 Seurat_setup.R](https://github.com/nyuhuyang/scRNAseq-MouseSkin/blob/master/R/Seurat_setup.R)

This script uses Seurat 2 to read the result from the Cell ranger outputs, perform normalization, scaling, dimension reduction, and unsupervised-clustering. The output Seurat object will be saved in the file `data/MouseSkin_{date}.Rda`

#### **2-3. Identify cell types**

[2 Identify_Cell_Types_Manually.R](https://github.com/nyuhuyang/scRNAseq-MouseSkin/blob/master/R/Identify_Cell_Types_Manually.R)

This script uses predefined cell type markers to identify cell types manually.
It contains source code for generating Fig 1B.

Besides, an R shiny app is built to identify cell types. [scRNAseq-MouseSkin](https://weillcornellmed.shinyapps.io/MouseSkin/)


[3 SingleR.R](https://github.com/nyuhuyang/scRNAseq-MouseSkin/blob/master/R/SingleR.R)

This script uses the SingleR package to identify cell types based on reference datasets.

#### **4-6. Generate Figures and data exploration**

[4 Major_Figures.R](https://github.com/nyuhuyang/scRNAseq-MouseSkin/blob/master/R/Major_Figures.R)

This script contains source code for generating Fig 1A, 1C, and 1S.

[5 Differential_analysis.R](https://github.com/nyuhuyang/scRNAseq-MouseSkin/blob/master/R/Differential_analysis.R), conducting differential analysis between Young and Aged sample.

[6 Monocle.R](https://github.com/nyuhuyang/scRNAseq-Immunosurveillance/blob/master/R/FGESA.R), performing Monocle analysis. 