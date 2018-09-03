#Identification of important genes
list.of.packages <- c("devtools","dplyr","pheatmap","VGAM", "irlba",
                      "matrixStats", "igraph", "combinat", "fastICA",
                      "grid", "reshape2", "plyr", "parallel", "methods")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
#check package
library("devtools")
library("DDRTree")
library("pheatmap")
library("M3Drop")
library("monocle")
library("reshape")
source("../R/Seurat_functions.R")
# 5.1 Importing data from Seurat object=================
lnames = load(file = "./output/MouseSkin_20180903.Rda");lnames
MouseSkin_Mono <- importCDS(MouseSkin, import_all = TRUE)

# 5.1.1 Estimate size factors and dispersions
# estimateSizeFactors() and estimateDispersions() will only work,
# and are only needed, if you are working with a CellDataSet 
# with a negbinomial() or negbinomial.size() expression family.
MouseSkin_Mono <- estimateSizeFactors(MouseSkin_Mono)
MouseSkin_Mono <- estimateDispersions(MouseSkin_Mono)

# 5.1.2 Filtering low-quality cells
MouseSkin_Mono <- detectGenes(MouseSkin_Mono, min_expr = 0.1)
print(head(fData(MouseSkin_Mono)))
print(head(pData(MouseSkin_Mono)))

# 5.1.3 If you are using RPC values to measure expresion, 
# as we are in this vignette, it's also good to look at the distribution
# of mRNA totals across the cells:
pData(MouseSkin_Mono)$Total_mRNAs <- Matrix::colSums(exprs(MouseSkin_Mono))
upper_bound <- 10^(mean(log10(pData(MouseSkin_Mono)$Total_mRNAs)) + 
                           2*sd(log10(pData(MouseSkin_Mono)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(MouseSkin_Mono)$Total_mRNAs)) - 
                           2*sd(log10(pData(MouseSkin_Mono)$Total_mRNAs)))
qplot(Total_mRNAs, data = pData(MouseSkin_Mono), color = conditions, geom = "density") +
        geom_vline(xintercept = lower_bound) +
        geom_vline(xintercept = upper_bound)

# 5.2 Classifying and counting cells of different types
# 5.2.1 Classifying cells with CellTypeHierarchy
all(row.names(pData(MouseSkin_Mono)) == names(MouseSkin@ident))
pData(MouseSkin_Mono)$CellType <- MouseSkin@ident
table(pData(MouseSkin_Mono)$CellType)
pie <- ggplot(pData(MouseSkin_Mono), aes(x = factor(1), fill = factor(CellType))) +
        geom_bar(width = 1)
pie + coord_polar(theta = "y") +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# 5.3 Constructing Single Cell Trajectories
# 5.3.1 choosing genes that define progress
expressed_genes <- row.names(subset(fData(MouseSkin_Mono), num_cells_expressed >= 10))
length(expressed_genes)
diff_test_res <- differentialGeneTest(MouseSkin_Mono[expressed_genes,],
                                      fullModelFormulaStr = "~ CellType",
                                      cores = 4) #takes long time
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
MouseSkin_Mono <- setOrderingFilter(MouseSkin_Mono, ordering_genes)
plot_ordering_genes(MouseSkin_Mono)

# 5.3.2 reduce data dimensionality
#Now we're ready to try clustering the cells:.
plot_pc_variance_explained(MouseSkin_Mono, return_all = F) # norm_method = 'log',
MouseSkin_Mono <- reduceDimension(MouseSkin_Mono, max_components = 2,
                            method = 'DDRTree')
#Trajectory step 3: order cells along the trajectory
MouseSkin_Mono <- orderCells(MouseSkin_Mono)

g1 <- plot_cell_trajectory(MouseSkin_Mono, color_by = "CellType",cell_size = 3)
g2 <- plot_cell_trajectory(MouseSkin_Mono, color_by = "State",cell_size = 3)
cowplot::plot_grid(g1,g2)
