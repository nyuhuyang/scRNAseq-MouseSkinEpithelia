# **The aging skin microenvironment dictates stem cell behavior**

This project provides the code developed for the study in [Yejing Ge et al, The aging skin microenvironment dictates stem cell behavior. PNAS (2020)](https://www.pnas.org/content/early/2020/02/18/1901720117). It focuses on comparing Young vs. Aged mouse skin epithelium with single-cell RNA-seq.

## **Requirements**

* R (tested in R version 3.5.2 (2018-12-20) -- "Eggshell Igloo")
* R librarys: Seurat (v2.3), SingleR (github version)

## **Data**

Raw counts data and processed Seurat object will be released after the final publication.

## **Reproduce results**

#### **1. Data preprocess**
[1 Seurat_setup.R](https://github.com/nyuhuyang/scRNAseq-MouseSkinEpithelia/blob/master/R/Seurat_setup.R)

Move Cell ranger output foder under `data`.
This script uses Seurat 2 to read the result from the Cell ranger outputs, perform normalization, scaling, dimension reduction, and unsupervised-clustering. The output Seurat object will be saved in the file `data/MouseSkin_{date}.Rda`

#### **2-3. Identify cell types**

[2 Identify_Cell_Types_Manually.R](https://github.com/nyuhuyang/scRNAseq-MouseSkinEpithelia/blob/master/R/Identify_Cell_Types_Manually.R)

This script uses predefined cell type markers to identify cell types manually.
![](https://github.com/nyuhuyang/scRNAseq-MouseSkinEpithelia/blob/master/Figs/F1B_dotplot.jpeg)

Besides, an R shiny app is built to identify cell types. [scRNAseq-MouseSkinEpithelia](https://weillcornellmed.shinyapps.io/MouseSkinEpithelia/)

![](https://github.com/nyuhuyang/scRNAseq-MouseSkinEpithelia/blob/master/Figs/Rshiny.png)

[3 SingleR.R](https://github.com/nyuhuyang/scRNAseq-MouseSkinEpithelia/blob/master/R/SingleR.R)

This script uses the SingleR package to identify cell types based on reference datasets.

#### **4-6. Generate Figures and data exploration**

[4 Major_Figures.R](https://github.com/nyuhuyang/scRNAseq-MouseSkinEpithelia/blob/master/R/Major_Figures.R)

This script contains source code for generating Fig 1A, 1C, and 1S.

![](https://github.com/nyuhuyang/scRNAseq-MouseSkinEpithelia/blob/master/Figs/F1A_title.jpeg)
![](https://github.com/nyuhuyang/scRNAseq-MouseSkinEpithelia/blob/master/Figs/F1C_Epd%20marker%20genes.jpeg)
![](https://github.com/nyuhuyang/scRNAseq-MouseSkinEpithelia/blob/master/Figs/F1C_HF%20marker%20genes.jpeg)
![](https://github.com/nyuhuyang/scRNAseq-MouseSkinEpithelia/blob/master/Figs/F1SA.jpeg)

[5 Differential_analysis.R](https://github.com/nyuhuyang/scRNAseq-MouseSkinEpithelia/blob/master/R/Differential_analysis.R), conducting differential analysis between Young and Aged sample.

[6 Monocle.R](https://github.com/nyuhuyang/scRNAseq-Immunosurveillance/blob/master/R/FGESA.R), performing Monocle analysis. 
Visit [Wiki](https://github.com/nyuhuyang/scRNAseq-MouseSkinEpithelia/wiki/The-aging-skin-microenvironment-dictates-stem-cell-behavior) for details.
