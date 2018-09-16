########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(SingleR)
library(dplyr)
library(ggrepel)
source("./R/Seurat_functions.R")
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here

MouseSkin_dgCMatrix <- Read10X(data.dir = paste0("./data/skin1/outs/raw_gene_bc_matrices/mm10/"))
MouseSkin_raw <- list("young" = MouseSkin_dgCMatrix[,grepl("-1",colnames(MouseSkin_dgCMatrix))],
                      "old" = MouseSkin_dgCMatrix[,grepl("-2",colnames(MouseSkin_dgCMatrix))])
MouseSkin_Seurat <- list()
samples <- c("young","old")
conditions <- c("young","old")
# 1 less, 2 more
# 1 yound, 2 old
for(i in 1:length(samples)){
    colnames(MouseSkin_raw[[i]]) <- paste0(samples[i],"_",colnames(MouseSkin_raw[[i]]))
    MouseSkin_Seurat[[i]] <- CreateSeuratObject(MouseSkin_raw[[i]],
                                           min.cells = 3,
                                           min.genes = 200,
                                           names.delim = "_")
    MouseSkin_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
MouseSkin_Seurat <- lapply(MouseSkin_Seurat, FilterCells, 
                            subset.names = "nGene", 
                            low.thresholds = 200, 
                            high.thresholds = Inf) %>% 
        lapply(NormalizeData) %>% 
        lapply(ScaleData) %>% 
        lapply(FindVariableGenes, do.plot = FALSE)

# we will take the union of the top 1k variable genes in each dataset for
# alignment note that we use 1k genes in the manuscript examples, you can
# try this here with negligible changes to the overall results
genes.use <- lapply(MouseSkin_Seurat, function(x) head(rownames(x@hvg.info), 1000))
genes.use <- unique(unlist(genes.use))
for(i in 1:length(samples)){
        genes.use <- intersect(genes.use, rownames(MouseSkin_Seurat[[i]]@scale.data))
}
length(genes.use) # 1/10 of total sample size


#======1.2 Perform a canonical correlation analysis (CCA) =========================
# run a canonical correlation analysis to identify common sources
# of variation between the two datasets.

MouseSkin <- RunCCA(MouseSkin_Seurat[[1]], MouseSkin_Seurat[[2]], 
                    genes.use = genes.use,num.cc = 50)
save(MouseSkin, file = "./data/MouseSkin_alignment.Rda")
remove(MouseSkin_raw,MouseSkin_Seurat)
GC()
# Calculate median UMI per cell
MouseSkin_raw_data <- as.matrix(x = MouseSkin@raw.data)
MouseSkin_young <- MouseSkin_raw_data[,grepl("young_",colnames(MouseSkin_raw_data))]
MouseSkin_old <- MouseSkin_raw_data[,grepl("old_",colnames(MouseSkin_raw_data))]
mean(colSums(MouseSkin_young)); mean(colSums(MouseSkin_old))
median(colSums(MouseSkin_young)); median(colSums(MouseSkin_old))
par(mfrow = c(1,2))
boxplot(colSums(MouseSkin_young), ylim = c(0, 21000))
title(main = "nUMI per cell in young mouse skin sample")
boxplot(colSums(MouseSkin_old),ylim = c(0, 21000))
title(main = "nUMI per cell in old mouse skin sample")
remove(MouseSkin_raw_data,MouseSkin_young,MouseSkin_old)
GC()
# CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = MouseSkin, reduction.use = "cca", group.by = "conditions", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = MouseSkin, features.plot = "CC1", group.by = "conditions", 
              do.return = TRUE)
plot_grid(p1, p2)

PrintDim(object = MouseSkin, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

p3 <- MetageneBicorPlot(MouseSkin, grouping.var = "conditions", dims.eval = 1:50, 
                        display.progress = TRUE, smooth = TRUE) # run on cluster
p3 + geom_smooth(method = 'loess')

DimHeatmap(object = MouseSkin, reduction.type = "cca", cells.use = 500, dim.use = c(1:3,25:20), 
           do.balanced = TRUE)

#======1.3 QC =========================
MouseSkin <- CalcVarExpRatio(object = MouseSkin, reduction.type = "pca",
                        grouping.var = "orig.ident", dims.use = 1:20)
MouseSkin
MouseSkin <- SubsetData(MouseSkin, subset.name = "var.ratio.pca",accept.low = 0.5) #10291 out of 10350
MouseSkin

mito.genes <- grep(pattern = "^mt-", x = rownames(x = MouseSkin@data), value = TRUE)
percent.mito <- Matrix::colSums(MouseSkin@raw.data[mito.genes, ])/Matrix::colSums(MouseSkin@raw.data)
MouseSkin <- AddMetaData(object = MouseSkin, metadata = percent.mito, col.name = "percent.mito")
#MouseSkin <- ScaleData(object = MouseSkin, genes.use = genes.use, display.progress = FALSE, 
#                         vars.to.regress = "percent.mito")
#Now we can run a single integrated analysis on all cells!
VlnPlot(object = MouseSkin, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 1)
par(mfrow = c(2, 1))
GenePlot(object = MouseSkin, gene1 = "nUMI", gene2 = "percent.mito", ylim = c(0, 0.6))
GenePlot(object = MouseSkin, gene1 = "nUMI", gene2 = "nGene", ylim = c(0, 4000))

MouseSkin <- FilterCells(object = MouseSkin, subset.names = c("nGene", "percent.mito"), 
                         low.thresholds = c(500, -Inf), high.thresholds = c(3500, 0.20))
GenePlot(object = MouseSkin, gene1 = "nUMI", gene2 = "percent.mito", 
         xlim = c(0,20000), ylim = c(0, 0.6))
GenePlot(object = MouseSkin, gene1 = "nUMI", gene2 = "nGene", 
         xlim = c(0,20000), ylim = c(0, 4000))
#======1.4 align seurat objects =========================
#Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned
set.seed(42)
MouseSkin <- AlignSubspace(object = MouseSkin, reduction.type = "cca", grouping.var = "conditions", 
                      dims.align = 1:20)
MouseSkin@project.name <- "Fuchs"
#Now we can run a single integrated analysis on all cells!

MouseSkin <- FindClusters(object = MouseSkin, reduction.type = "cca.aligned", dims.use = 1:20, 
                     resolution = 1.0, force.recalc = T, save.SNN = TRUE)

MouseSkin <- RunTSNE(object = MouseSkin, reduction.use = "cca.aligned", dims.use = 1:20, 
                do.fast = TRUE)

p1 <- TSNEPlot(MouseSkin, do.return = T, pt.size = 1, group.by = "orig.ident")
p2 <- TSNEPlot(MouseSkin, do.return = T, pt.size = 1, group.by = "ident")
#png('./output/TSNESplot_alignment.png')
plot_grid(p1, p2)
#dev.off()
TSNEPlot(object = MouseSkin,do.label = T, group.by = "ident", 
         do.return = TRUE, no.legend = T,
         pt.size = 1,label.size = 8 )+
        ggtitle("TSNE plot of young and old mouse skin")+
        theme(text = element_text(size=20),     							
              plot.title = element_text(hjust = 0.5))

save(MouseSkin, file = "./data/MouseSkin_alignment.Rda")
