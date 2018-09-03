library(Seurat)
library(MAST)
library(dplyr)
library(tidyr)
library(kableExtra)
source("../R/Seurat_functions.R")

#====== 2.1 identify phenotype for each cluster  ==========================================
lnames = load(file = "./data/MouseSkin_alignment.Rda")
lnames
Featureplot <- function(x,object = MouseSkin,...){
    p <- FeaturePlot(object = object, 
                     reduction.use = "tsne",
                     features.plot = x, min.cutoff = NA, 
                     cols.use = c("lightgrey","blue"), pt.size = 1,...)
    return(p)
}

Reference_cell_populations_sig_genes <- readxl::read_excel("./doc/Reference_cell_populations_sig_genes_HY.xlsx")
Skin_Markers <- as.data.frame(Reference_cell_populations_sig_genes[,-1])
colnames(Skin_Markers) = Skin_Markers[1,]
Skin_Markers = Skin_Markers[-1,]
Skin_Markers = apply(Skin_Markers, 2, function(x) sub(" \\(.*","",x))
Skin_Markers = as.data.frame(Skin_Markers[1:15,])
Skin_Markers = df2list(Skin_Markers) # must convert to list to exclude NA
Skin_Markers = lapply(Skin_Markers,function(x) MouseGenes(MouseSkin,x)) # must genes exist in MouseSkin
Featureplot(Skin_Markers$Proliferative)
#  test Find all ConservedMarkers=============
FindAllConservedMarkers <- function(object, grouping.var, assay.type = "RNA", 
                                    meta.method = minimump, print.bar = FALSE){
        idents <- as.data.frame(table(object@ident))
        idents = as.character(idents$Var1)
        Conserved.Markers <- list()
        for(i in 1:length(idents)){
                Conserved.Markers[[i]] <- FindConservedMarkers(object,ident.1 = idents[i], 
                                                                grouping.var = grouping.var, 
                                                                print.bar = print.bar)
                Conserved.Markers[[i]]$cluster <- idents[i]
                Conserved.Markers[[i]]$gene <- rownames(Conserved.Markers[[i]])
                rownames(Conserved.Markers[[i]]) = NULL
        }
        all.Conserved.Markers <- do.call("rbind", Conserved.Markers)
        return(all.Conserved.Markers)
}

all.Conserved.Markers <- FindAllConservedMarkers(MouseSkin,grouping.var = "conditions")
write.csv(all.Conserved.Markers,"./output/Conserved_Markers_by_cluster.csv")

#AllMarkers_cluster <- FindAllMarkers.UMI(MouseSkin,test.use = "MAST")
#write.csv(AllMarkers_cluster,"./output/AllMarkers_by_cluster.csv")

df_Skin_Markers <- list2df(Skin_Markers)
df_Skin_Markers <- gather(df_Skin_Markers, key = "cell_type", value = "gene",
               na.rm = T)
AllMarkers.cluster <- inner_join(all.Conserved.Markers, df_Skin_Markers,by = "gene")
count.cluster <- AllMarkers.cluster[,c("cluster","cell_type")]
count.cluster = count.cluster %>% group_by(cluster) %>% dplyr::count(cell_type)
write.csv(count.cluster,"./output/Count_by_cluster.csv")
count.cluster = count.cluster %>% group_by(cluster) %>% top_n(2, n)


#====== 2.2 Rename ident ==========================================
# 
lnames = load(file = "./data/MouseSkin_alignment.Rda")
lnames
table(MouseSkin@ident)
idents <- as.data.frame(table(MouseSkin@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Basal Epd",
                     "Basal HF",
                     "Supra-basal Epd",
                     "Basal Epd",
                     "Basal HF",
                     "Bulge (HF-SCs)",
                     "Proliferative",
                     "Hair Germ",
                     "Basal Epd",
                     "Supra-basal HF",
                     "Basal Epd",
                     "Sebaceous Gland",
                     "Basal HF")

MouseSkin@ident <- plyr::mapvalues(x = MouseSkin@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
# TSNEPlot
MouseSkin@ident <- factor(MouseSkin@ident, levels = (c("Hair Germ",
                                                       "Basal HF",
                                                       "Proliferative",
                                                       "Basal Epd",
                                                       "Supra-basal Epd",
                                                       "Sebaceous Gland",
                                                       "Bulge (HF-SCs)",
                                                       "Supra-basal HF")))
TSNEPlot.1(object = MouseSkin,do.label = T, group.by = "ident", 
           do.return = TRUE, no.legend = T, 
           colors.use = c(singler.colors[c(1,3:8)],"#33b3a6"),
           pt.size = 1,label.size = 5,label.repel = T,force=1)+
        ggtitle("TSNEPlot of all clusters")+
        theme(text = element_text(size=20),							
              plot.title = element_text(hjust = 0.5,size=18, face = "bold"))

# dotplot
MouseSkin@ident <- factor(MouseSkin@ident, levels = (c("Hair Germ",
                                                       "Bulge (HF-SCs)",
                                                       "Basal HF",
                                                       "Supra-basal HF",
                                                       "Sebaceous Gland",
                                                       "Proliferative",
                                                       "Basal Epd",
                                                       "Supra-basal Epd")))
markers.to.plot = unique(c(Skin_Markers$`Suprabasal Epidermis (SB Epd)`,
                    Skin_Markers$`Basal Epidermis (Basal Epd)`,
                    Skin_Markers$Proliferative,
                    Skin_Markers$`Sebaceous Gland`,
                    Skin_Markers$`Suprabasal Hair Follicle (SB HF)`,
                    Skin_Markers$`Basal Hair Follicle (Basal HF)`,
                    Skin_Markers$`Bulge (Hair Follicle Stem Cells)`,
                    Skin_Markers$`Hair Germ (Hair Follicle Progenitors)`))

# dot Plots 
markers.to.plot <- c("Mt4","Krt1","Krt77","Krt10","Krtdap","Calm4",
                     "Gata3","Ly6a","Klf5","Ndufa4l2","Serpinb2","Pou3f1","Ifngr1","Il1r2","Il6ra","Il33",
                     Skin_Markers$Proliferative,
                     Skin_Markers$`Sebaceous Gland`[-2],
                     Skin_Markers$`Suprabasal Hair Follicle (SB HF)`[-1],
                     "Krt17","Aqp3","Aldh3a1","Sostdc1","Gstm1","Gstm5","Cd200","Fst","Calml3","Sox9","Ctgf",
                     Skin_Markers$`Bulge (Hair Follicle Stem Cells)`[-c(5:7,10,11)],
                     "Tgm5","Lhx2","Lgr5","Cldn10")

sdp <- SplitDotPlotGG(MouseSkin, genes.plot = rev(markers.to.plot), 
                      cols.use = c("blue","red"), x.lab.rot = T,plot.legend = T,
                      dot.scale = 8, do.return = T, grouping.var = "conditions")

#====== 2.3 Cell number and % ==========================================
table(MouseSkin@ident, MouseSkin@meta.data$orig.ident) %>% kable() %>% kable_styling()
table(MouseSkin@ident, MouseSkin@meta.data$orig.ident) %>% 
        prop.table(margin = 2) %>% kable() %>% kable_styling()# margin = 2: divided by the sum of the column cells

#====== 2.4 Vln Plot==========================================

VlnPlot(object = pbmc, features.plot = c("MS4A1", "CD79A"))

Reference_cell_populations_sig_genes <- read_excel("doc/Reference cell populations sig genes.xlsx", 
                                                   +     sheet = "Violin plot genes")