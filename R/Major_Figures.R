library(Seurat)
library(dplyr)
library(gridExtra)
library(magrittr)
library(ggplot2)
library(cowplot)
library(fgsea)
library(tibble)
library(ggsci)
source("R/Seurat3_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

(load(file = "data/MouseSkin_20180903.Rda"))
object <- UpdateSeuratObject(MouseSkin)
remove(MouseSkin);GC()
object[["major_cell.types"]] = Idents(object)


# rerun FindClusters
set.seed(101)
object %<>% FindNeighbors(reduction = "cca.aligned",dims = 1:20)
object %<>% FindClusters(resolution = 1)
object %<>% FindClusters(resolution = 1.2)
Idents(object) ="RNA_snn_res.1.2"
object %<>% RenameIdents('0' = 'Basal Epd',
                         '1' = 'Infund',
                         '2' = 'SupraBasal Epd',
                         '3' = 'Basal Epd',
                         '4' = 'Isthmus',
                         '5' = 'Bulge',
                         '6' = 'Proliferation',
                         '7' = 'Hair germ',
                         '8' = 'Basal Epd',
                         '9' = 'SupraBasal HF',
                         '10'= 'Basal Epd',
                         '11'= 'SebGland',
                         '12'= 'Bulge',
                         '13'= 'Basal Epd')
cell.types = c("Bulge","Hair germ","Isthmus","Infund",
               "SebGland","SupraBasal HF","Basal Epd",
               "SupraBasal Epd","Proliferation")
levels(object) = cell.types
object@meta.data[,"cell.types"] = Idents(object)
if(!"tSNE_1" %in% colnames(object@meta.data)) object@meta.data = cbind(object@meta.data, 
                                                        object@reductions$tsne@cell.embeddings)
object@meta.data[object$tSNE_2 >-19 & object$RNA_snn_res.1.2 %in% 8,"cell.types"] = 'Basal Epd'
Idents(object) = "cell.types"
# Figure S1B
TSNEPlot.1(object,label = F,no.legend = F,
           cols = rev(gg_color_hue(length(cell.types))),
           pt.size = 1.5,legend.size = 15, do.print = T)

# Figure 1A
object@meta.data[,"conditions"] %<>% gsub("old","Aged",.)
object@meta.data[,"conditions"] %<>% gsub("young","Young",.)
Idents(object) = "conditions"

levels(object) = c("Young","Aged")
object[["conditions"]] = Idents(object)
TSNEPlot.1(object,label = F,no.legend = F,group.by = "conditions",
           cols = rev(gg_color_hue(2)),
           pt.size = 1.5,do.print = T)

young <- subset(object, idents = "Young")
old <- subset(object, idents = "Aged")

NoAxisTitle <- theme(axis.text.x = element_text(size=8),
                     axis.text.y = element_blank(),
                     axis.title.x = element_blank(), axis.title.y = element_blank())
g1 <- TSNEPlot.1(young,label = F,no.legend = T,group.by = "cell.types",
                 cols = rev(gg_color_hue(length(cell.types))),
                 label.size = 5,
                 pt.size = 0.5,do.print = F) +NoAxisTitle+
        xlim(-50,45)
g2 <- TSNEPlot.1(old,label = F,no.legend = T,group.by = "cell.types",
                 cols = rev(gg_color_hue(length(cell.types))),
                 pt.size = 0.5,do.print = F)+ NoAxisTitle
g3 <- TSNEPlot.1(object,label = F,no.legend = T,group.by = "conditions",
                 cols = rev(gg_color_hue(2)),
                 pt.size = 0.5,do.print = F)+ NoAxisTitle
g4 <- TSNEPlot.1(object,label = F,no.legend = T,group.by = "cell.types",
                 cols = rev(gg_color_hue(length(cell.types))),
                 pt.size = 0.5,do.print = F)+ NoAxisTitle
jpeg(paste0(path,"F1A.jpeg"),units="in", width=10, height=2.5,res=600)
grid.arrange(g1, g2, g3, g4, nrow = 1)
dev.off()

g1 <- TSNEPlot.1(young,label = F,no.legend = T,group.by = "cell.types",
                 cols = rev(gg_color_hue(length(cell.types))),
                 label.size = 5,
                 title= "Young",
                 pt.size = 0.5,do.print = F) +
        NoAxisTitle+
        xlim(-50,45)
g2 <- TSNEPlot.1(old,label = F,no.legend = T,group.by = "cell.types",
                 cols = rev(gg_color_hue(length(cell.types))),
                 title= "Aged",
                 pt.size = 0.5,do.print = F)+ 
        NoAxisTitle
g3 <- TSNEPlot.1(object,label = F,no.legend = T,group.by = "conditions",
                 cols = rev(gg_color_hue(2)),
                 title= "Young vs. Aged",
                 pt.size = 0.5,do.print = F)+ 
        NoAxisTitle
g4 <- TSNEPlot.1(object,label = F,no.legend = T,group.by = "cell.types",
                 cols = rev(gg_color_hue(length(cell.types))),
                 title= "Combine",
                 pt.size = 0.5,do.print = F)+ 
        NoAxisTitle
jpeg(paste0(path,"F1A_title.jpeg"),units="in", width=10, height=2.5,res=600)
grid.arrange(g1, g2, g3, g4, nrow = 1)
dev.off()


# 4) Violin plot for selected HF or Epd genes (2nd sheet on the excel file) 
VioPlot_genes <- c("Krt24","Sox9","Lgr5","Grem1","Lhx2","Tnc",
                   "Ly6a","Klf5","Gata3","Tfap2c","Krt10","Grhl1")
Vioplot_genes = FilterGenes(object,VioPlot_genes)
Idents(object) = "cell.types"
Bulge <- subset(object, idents = 'Bulge')
Idents(Bulge) = "conditions"

jpeg(paste0(path,"F1C_HF marker genes.jpeg"),units="in", width=10, height=7,res=600)
VlnPlot(object = Bulge, features = Vioplot_genes[1:6],group.by = "conditions")
dev.off()

jpeg(paste0(path,"F1C_Epd marker genes.jpeg"),units="in", width=10, height=7,res=600)
VlnPlot(object = Bulge, features = Vioplot_genes[7:12],group.by = "conditions")
dev.off()