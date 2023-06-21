#clear objects in workspace
rm(list = ls())
options(stringsAsFactors = F)

#required packages
library("Seurat")
library("dplyr")
library("ggplot2")
library(RColorBrewer)
library(cowplot)
library(reticulate)
library(Matrix)

#store rds file
saveRDS(B_all, "/singleron/april29_matrix/B_newest_analyzed.rds")

#read stored rds files
B_all <- readRDS("/publication/analyzed/B.rds")

#gene of interest expression in sample
Idents(B_all) <- "cellcluster"
jpeg("/singleron/test_B.jpg", width = 2400, height = 2400)
FeaturePlot(B_all, reduction = "umap", pt.size = 5, features = c("bf3-583"))
dev.off()

#run integrated analysis of all rds files
all <- merge(B_all,C_all)
all <- merge(all, G_all)
all <- merge(all, D_all)
all <- merge(all, N2_all)
all <- merge(all, E_all)
all <- merge(all, F_all)
table(all@meta.data$orig.ident)
all.list <- SplitObject(all, split.by = "orig.ident")
all.list <- all.list[c("B1","B2","C_1","C_2","C_3","G_sample1","G_sample2",
                     "D_sample4","D_sample1","D_sample2","D_sample3",
                     "BW-1","BW-2","E_sample1","E_sample2","E_sample3",
                     "E_sample4","F_sample1","F_sample2","F_sample3",
                     "F_sample4","F_sample5")]
all.list <- lapply(X = all.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = all.list)
all.list <- lapply(X = all.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
all.anchors <- 
  FindIntegrationAnchors(object.list = all.list, anchor.features = 2500, 
                         reduction = "rpca", k.filter = 200, dims = 1:50)
all.combined <- IntegrateData(anchorset = all.anchors)
DefaultAssay(all.combined) <- "integrated"
all.genes <- rownames(all.combined)
all.combined <- ScaleData(all.combined, features = all.genes, verbose = FALSE)
all.combined <- 
  RunPCA(all.combined, features = VariableFeatures(object = all.combined), 
         verbose = FALSE)
all.combined <- RunUMAP(all.combined, dims = 1:50)
DefaultAssay(all.combined) <- "RNA"