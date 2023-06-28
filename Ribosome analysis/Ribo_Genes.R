library(Seurat)
library(sctransform)
library(dplyr)
library(tidyverse)
library(readxl)
library(stringr)

GerberAxo <- readRDS(file = "NewAnnotation/Gerber_Seurat_NewAnno_GeneID_prrx1.rds")

DimPlot(GerberAxo, reduction = "umap", 
        pt.size = 0.1, label = T, label.size = 6, group.by = "Stage") + NoLegend()

Ribo.genes <- read.csv(file = "Ribosome biogenesis genes_ID_Final.csv", header = T)


levels(GerberAxo@meta.data$Stage)[levels(GerberAxo@meta.data$Stage)=='Uninjured_1apa'] <- 'Uninjured_1dpa'


GerberAxo@meta.data$Stage <- factor(GerberAxo@meta.data$Stage, 
                                    levels = c("Uninjured_0dpa", 
                                               "Uninjured_1dpa", 
                                               "Uninjured_3mpa",
                                               "Blastema_3dpa",
                                               "Blastema_5dpa",
                                               "Blastema_8dpa",
                                               "Blastema_11dpa",
                                               "Blastema_18dpa"
                                               ))



GerberAxo <- AddModuleScore(object = GerberAxo, 
                            features = list(Ribo.genes$GeneID), 
                            name = "Ribo.score")
head(GerberAxo@meta.data)
FeaturePlot(GerberAxo, features = "Ribo.score1")
VlnPlot(GerberAxo, features = "Ribo.score1", group.by = "Stage")
