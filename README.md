# **The aging skin microenvironment dictates stem cell behavior**

This project provides the code developed in the study of Yejing Ge _et al._ **_"The aging skin microenvironment dictates stem cell behavior"_**. It focuses on comparing Young vs. Aged mouse skin with single-cell RNA-seq.

## **Requirements**

* R (tested in R version 3.5.2 (2018-12-20) -- "Eggshell Igloo")
* R librarys: Seurat (v2.3), SingleR (github version)

## **Data**

Raw counts data and processed Seurat object will be released after final publication.

## **Reproduce results**

#### **1. Data preprocess**
Move Cell ranger output foder under `data`.
[1 Seurat_setup.R](https://github.com/nyuhuyang/scRNAseq-MouseSkin/blob/master/R/Seurat_setup.R)

This script uses Seurat 2 to read the result from the Cell ranger outputs; perform normalization, scaling, dimension reduction, and unsupervised-clustering. The output Seurat object will be saved in the file `data/MouseSkin_{date}.Rda`

#### **2-3. Identify cell types**

[2 Identify_Cell_Types_Manually.R](https://github.com/nyuhuyang/scRNAseq-MouseSkin/blob/master/R/Identify_Cell_Types_Manually.R)

This script use predefinde cell type markers to manually identify cell types.

In addition, an R shiny app is built to identify cell types. [scRNAseq-MouseSkin]()


[3 SingleR.R](https://github.com/nyuhuyang/scRNAseq-MouseSkin/blob/master/R/SingleR.R)

This script use SingleR package to identify cell types based reference datasets.

#### **6-9. Characterize CD45+ tumor infiltrate of mouse tumors**

[6 Figures.R](https://github.com/nyuhuyang/scRNAseq-Immunosurveillance/blob/master/R/Figures.R), removing unwanted cells, relabel the clustering, generate tsne plot and bar plots. 

For more information see
[_Generate TSNE plots and compare gene expression in individual cell types (figure 8a-b)_](https://github.com/nyuhuyang/scRNAseq-Immunosurveillance/wiki/1.-Generate-TSNE-plots-and-compare-gene-expression-in-individual-cell-types)

[7 Differential_analysis.R](https://github.com/nyuhuyang/scRNAseq-Immunosurveillance/blob/master/R/Differential_analysis.R), conducting differential analysis between control and NAM treated tumor sample, and generate heatmaps.

[8 FGESA.R](https://github.com/nyuhuyang/scRNAseq-Immunosurveillance/blob/master/R/FGESA.R), performing gene set enrichment analysis based on the results from previous step 7. 

For more information, see [_Differential analysis and gene set enrichment analysis (figure S7, S8)_](https://github.com/nyuhuyang/scRNAseq-Immunosurveillance/wiki/2.-Differential-analysis-and-gene-set-enrichment-analysis)
