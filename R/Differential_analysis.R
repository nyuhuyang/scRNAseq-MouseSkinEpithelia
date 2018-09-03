########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(kableExtra)
library(SingleR)
library(MAST)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")

#4.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 4.1.1 load data=====
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

# 4.2  for 2018/08/31 email ==================
Reference_cell_populations_sig_genes <- readxl::read_excel("./doc/Reference_cell_populations_sig_genes_HY.xlsx")
Skin_Markers <- as.data.frame(Reference_cell_populations_sig_genes[,-1])
colnames(Skin_Markers) = Skin_Markers[1,]
Skin_Markers = Skin_Markers[-1,]
Skin_Markers = apply(Skin_Markers, 2, function(x) sub(" \\(.*","",x))
Skin_Markers = as.data.frame(Skin_Markers[1:15,])
Skin_Markers = df2list(Skin_Markers) # must convert to list to exclude NA
Skin_Markers = lapply(Skin_Markers,function(x) MouseGenes(MouseSkin,x)) # must genes exist in MouseSkin
df_Skin_Markers <- list2df(Skin_Markers)
df_Skin_Markers <- gather(df_Skin_Markers, key = "cell_type", value = "gene",
                          na.rm = T)

# 4.2.1 VioPlot_genes=====
VioPlot_genes <- readxl::read_excel("./doc/Reference_cell_populations_sig_genes_HY.xlsx", 
               sheet = "Violin plot genes")
VioPlot_genes = as.vector(as.matrix(VioPlot_genes[,1]))
df_Skin_Markers[(df_Skin_Markers$gene %in% VioPlot_genes),"cell_type"]
Vioplot_genes = MouseGenes(MouseSkin,VioPlot_genes)
VlnPlot(object = MouseSkin, features.plot = Vioplot_genes[31:32],x.lab.rot = T)

# 4.2.2 DE analysis for a) young vs old Bulge,Basal Epidermis, Basal Hair Follicle
FindAllMarkersbyAge<- function(object, ident.use ,colors.use = NULL, 
                               select.plots = c(2,1),pt.size = 1,label.size = 5,
                               label.repel = T,force=1,
                               test.use = "MAST",...){
        object <- SubsetData(object, ident.use = ident.use)
        SplitTSNEPlot(object = object, do.label = T, 
                      colors.use = colors.use, select.plots = select.plots,
                      pt.size = pt.size,label.size = label.size,
                      label.repel = label.repel,force=force)
        ident.new = paste(ident.use, object@meta.data$orig.ident,sep = "_")
        object <- SetIdent(object, ident.use = ident.new)
        print(table(object@ident))
        gde.all <- FindAllMarkers.UMI(object,test.use = test.use,...)
        gde.all =  dplyr::select(gde.all, "gene",everything())
        rownames(gde.all) = NULL
        
        path <- paste("./output",gsub("-","",Sys.Date()),sep = "/")
        dir.create(path, recursive = T)
        write.csv(gde.all,file= paste0(path,"/",ident.use,"_old_vs_young.csv"))
        return(gde.all)
}

BulgebyAge = FindAllMarkersbyAge(MouseSkin,ident.use = "Bulge (HF-SCs)" ,
                                 colors.use = "#1B9E77")
BasalEpdbyAge = FindAllMarkersbyAge(MouseSkin,ident.use = "Basal Epd" ,
                                 colors.use = "#F0027F")
BasalHFbyAge = FindAllMarkersbyAge(MouseSkin,ident.use = "Basal HF" ,
                                    colors.use = "#FDC086")
# find and print differentially expressed genes across conditions ================
# combine SetIdent,indMarkers and avg_UMI

# 4.2.3 To Monocle
#source("http://bioconductor.org/biocLite.R")
#biocLite()
biocLite("monocle")                  
library(monocle)

MouseSkin_monocle = importCDS(MouseSkin)
