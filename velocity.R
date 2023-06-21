#required packages
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(reticulate)
library(anndata)
library(reshape)
library("RColorBrewer")
library("pheatmap")

#prepare analyzed bam file
E_3 <- subset(E_all, subset = orig.ident == "E_sample3")
E_3$velocity_Barcodes <- paste0("CB:Z:",E_3$Barcodes)
write.table(E_3$velocity_Barcodes, 
            "/singleron/velocity/E3_filtered_barcodes.txt", 
            append = FALSE, sep = "\t", dec = ".", 
            row.names = FALSE, col.names = FALSE)
samtools view -H E1-3_Aligned.sortedByCoord.out.bam.featureCounts.bam | wc -l
samtools view -h E1-3_Aligned.sortedByCoord.out.bam.featureCounts.bam | head -n 10105 | samtools view -bS - > test_E1_3.bam
samtools view -h test_E1_3.bam
samtools view -h cellsorted_test_E1_3.bam
velocyto run -o output \
test_E1_3.bam /reference/1109_updated_reference/bf.2111.new.gtf \
velocyto run -b test_filtered_barcodes.txt -o output \
test_E1_3.bam /reference/1109_updated_reference/bf.2111.new.gtf \
samtools view -H test_E1_3.bam > SAM_header
samtools view test_E1_3.bam | LC_ALL=C grep -F -f test_filtered_barcodes.txt > filtered_SAM_body
cat SAM_header filtered_SAM_body > filtered.sam
samtools view -b filtered.sam > filtered_test_E1_3.bam

#prepare loom file
velocyto run -b E3_test_barcodes.txt -o output \
E1-3_Aligned.sortedByCoord.out.bam.featureCounts.bam \
/reference/1109_updated_reference/bf.2111.new.gtf
#merge different samples/lanes in the same loom file in python
import loompy
files = ["E_combined.loom","F_combined.loom"]
loompy.combine(files, "EF_combined.loom", key="Accession")

#required packages in R
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(ggplot2)

#Read combined rds files for two developmental stages of interest
EF_all <- readRDS("/publication/N4_T1_integrated.rds")
efldat <- ReadVelocity(file = "/singleron/velocity/EF_combined.loom")
efm <- as.Seurat(x = efldat)
#Dataset contains two assays: the number of spliced and unspliced reads
efm@assays
#Generate the UMAP plot, but save it as an object
Idents(EF_all) <- "cellcluster"
jpeg("/publication/supp18_ef_full_all_umap.jpg", width = 3000, height = 3000)
DimPlot(EF_all, reduction = "umap", 
        pt.size = 1, label = TRUE, repel = TRUE, label.size = 10) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  NoLegend()
dev.off()
p <- UMAPPlot(EF_all) 
pbuild <- ggplot2::ggplot_build(p)
pdata <- pbuild$data[[1]] # Pull the data used for the plot
length(pdata$colour)
#Subset cells included in Seurat
EF_all$velocity_barcode <- paste(EF_all$orig.ident, EF_all$Barcodes, sep = "_")
EF_all$velocity_barcode <- gsub("E_sample1_","E1-1_Aligned_F4TJI:",
                                EF_all$velocity_barcode)
EF_all$velocity_barcode <- gsub("E_sample2_","E1-2_Aligned_VD8ET:",
                               EF_all$velocity_barcode)
EF_all$velocity_barcode <- gsub("E_sample3_","E_sample3_Aligned_1I10X:",
                               EF_all$velocity_barcode)
EF_all$velocity_barcode <- gsub("E_sample4_","E1-4_Aligned_KQY8D:",
                               EF_all$velocity_barcode)
EF_all$velocity_barcode <- gsub("F_sample1_","F1-1_Aligned_3R5D9:",
                                EF_all$velocity_barcode)
EF_all$velocity_barcode <- gsub("F_sample2_","F1-2_Aligned_Y32TF:",
                                EF_all$velocity_barcode)
EF_all$velocity_barcode <- gsub("F_sample3_","F2-1_Aligned_X82SE:",
                                EF_all$velocity_barcode)
EF_all$velocity_barcode <- gsub("F_sample4_","F2-2_Aligned_WV2KB:",
                                EF_all$velocity_barcode)
EF_all$velocity_barcode <- gsub("F_sample5_","F2-3_Aligned_MTINU:",
                                EF_all$velocity_barcode)
EF_all$velocity_barcode <- paste0(EF_all$velocity_barcode, "x")
subset_cells <- subset(efm, cells = EF_all$velocity_barcode)
subset_cells
subset_cells_barcode <- Cells(subset_cells)
sub_cells <- read.csv("/velocity/supp18_N4_T1_metadata.csv", header = T)
sub_cells_umap <- 
  read.csv("/velocity/supp18_N4_T1_cell_embeddings.csv", header = T)
sub_cells$umap1 <- sub_cells_umap$UMAP_1
sub_cells$umap2 <- sub_cells_umap$UMAP_2
sub_cells$color <- pdata$colour
sub_cells$velocity_barcode <- paste(sub_cells$orig.ident, sub_cells$Barcodes, 
                                    sep = "_")
sub_cells$velocity_barcode <- gsub("E_sample1_","E1-1_Aligned_F4TJI:",
                                   sub_cells$velocity_barcode)
sub_cells$velocity_barcode <- gsub("E_sample2_","E1-2_Aligned_VD8ET:",
                                   sub_cells$velocity_barcode)
sub_cells$velocity_barcode <- gsub("E_sample3_","E_sample3_Aligned_1I10X:",
                                   sub_cells$velocity_barcode)
sub_cells$velocity_barcode <- gsub("E_sample4_","E1-4_Aligned_KQY8D:",
                                   sub_cells$velocity_barcode)
sub_cells$velocity_barcode <- gsub("F_sample1_","F1-1_Aligned_3R5D9:",
                                   sub_cells$velocity_barcode)
sub_cells$velocity_barcode <- gsub("F_sample2_","F1-2_Aligned_Y32TF:",
                                   sub_cells$velocity_barcode)
sub_cells$velocity_barcode <- gsub("F_sample3_","F2-1_Aligned_X82SE:",
                                   sub_cells$velocity_barcode)
sub_cells$velocity_barcode <- gsub("F_sample4_","F2-2_Aligned_WV2KB:",
                                   sub_cells$velocity_barcode)
sub_cells$velocity_barcode <- gsub("F_sample5_","F2-3_Aligned_MTINU:",
                                   sub_cells$velocity_barcode)
sub_cells$velocity_barcode <- paste0(sub_cells$velocity_barcode, "x")
sub_cells_clean <- subset(sub_cells, 
                          (sub_cells$velocity_barcode %in% subset_cells_barcode))
barcodes <- sub_cells_clean$velocity_barcode
identity <- sub_cells_clean$cellcluster
class <- sub_cells_clean$celltype
group <- sub_cells_clean$cellgroup
color <- sub_cells_clean$color
umap1 <- sub_cells_clean$umap1
umap2 <- sub_cells_clean$umap2
Idents(object = subset_cells) <- "idents"
clusters <- read.csv("/publication/EF_metadata.csv", header = T)
newnames <- clusters$group
names(newnames) <- levels(subset_cells)
subset_cells <- RenameIdents(subset_cells, newnames)
subset_cells$group <- Idents(subset_cells)
subset_cells <- AddMetaData(object = subset_cells, 
                            metadata = identity, col.name = 'idents')
subset_cells <- AddMetaData(object = subset_cells,
                            metadata = class, col.name = 'class')
subset_cells <- AddMetaData(object = subset_cells, 
                            metadata = group, col.name = 'group')
subset_cells <- AddMetaData(object = subset_cells,
                            metadata = color, col.name = 'color')
subset_cells <- AddMetaData(object = subset_cells,
                            metadata = umap1, col.name = 'umap1')
subset_cells <- AddMetaData(object = subset_cells,
                            metadata = umap2, col.name = 'umap2')
Idents(object = subset_cells) <- "group"
sub_neural <- subset(subset_cells,
                     idents = c("Neuromesoderm","Epithelial_ectoderm",
                                "Neural_CNS","Neural_PNS"))
sub_mesoderm <- subset(subset_cells, idents = c("Mesoderm"))
Idents(object = subset_cells) <- "idents"
sub_single <- subset(subset_cells,idents = c("T1:61"))
sub_endoderm <- subset(subset_cells, idents = c("Endoderm"))
sub_meso <- merge(sub_mesoderm,sub_single)
Idents(object = sub_endoderm) <- "idents"
E_list <- unique(sub_endoderm$idents)
E_downsampled_list <- c("0")
for(name in E_list) {
  E_mectoder_name <- subset(sub_endoderm,idents = name)
  if (length(colnames(E_mectoder_name)) > 500){
    E_downsampled_name <- sample(colnames(E_mectoder_name), 
                                 size = 500, replace = F)
  }
  else{
    E_downsampled_name = colnames(E_mectoder_name)
  }
  E_downsampled_list <- c(E_downsampled_list, E_downsampled_name)
  print(paste("finished",name))
}
length(E_downsampled_list)
E_downsampled_list <- E_downsampled_list[-1]
cells_downsampled <- subset(subset_cells,cells = E_downsampled_list)
#Split the combined object into a list, with each dataset as an element
merged_list <- SplitObject(cells_downsampled, split.by = "orig.ident")
merged_list <- merged_list[c("E1-1", "E1-2", "E", "E1-4", 
                             "F1-1", "F1-2", "F2-1", "F2-2", "F2-3")]
merged_list <- lapply(X = merged_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = merged_list)
merged_list <- lapply(X = merged_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
merged_anchors <- 
  FindIntegrationAnchors(object.list = merged_list, anchor.features = 2500, 
                         reduction = "rpca", k.filter = 200, dims = 1:50)
#Return a Seurat object
merged_ect.integd <- IntegrateData(anchorset = merged_anchors, dims = 1:50)
DefaultAssay(merged_ect.integd) <- "integrated"
all.genes <- rownames(merged_ect.integd)
merged_ect.integd <- 
  ScaleData(merged_ect.integd, features = all.genes, verbose = FALSE)
merged_ect.integd <- 
  RunPCA(merged_ect.integd, 
         features = VariableFeatures(object = merged_ect.integd), 
         verbose = FALSE)
merged_ect.integd <- RunUMAP(merged_ect.integd, dims = 1:50)
Idents(object = merged_ect.integd) <- "idents"
merged_ect.integd[["umap"]]@cell.embeddings[,1] <- merged_ect.integd$umap1
merged_ect.integd[["umap"]]@cell.embeddings[,2] <- merged_ect.integd$umap2
jpeg("/publication/supp18_ef_endoderm_umap.jpg", width = 2000, height = 2000)
DimPlot(merged_ect.integd, reduction = "umap", 
        pt.size = 2, label = TRUE, repel = TRUE, label.size = 10) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  NoLegend()
dev.off()
DefaultAssay(merged_ect.integd) <- "ambiguous"
SaveH5Seurat(merged_ect.integd, filename = "EF_endoderm_500_0221.h5Seurat")
Convert("EF_endoderm_500_0221.h5Seurat", dest = "h5ad")
merged_ect.integd <- 
  RunVelocity(object = merged_ect.integd, deltaT = 1, kCells = 25, 
              fit.quantile = 0.02)
#Extra information
merged_ect.integd@tools$RunVelocity$projected[1:10, 1:10]
#Predicted gamma values for each gene
head(merged_ect.integd@tools$RunVelocity$gamma)
#kNN
merged_ect.integd@tools$RunVelocity$cellKNN[1:10, 1:10]
saveRDS(merged_ect.integd, "/velocity/EF_endoderm_500_0221.rds")
#Specify colors
ident.colors <- unique(merged_ect.integd$color)
names(x = ident.colors) <- levels(x = merged_ect.integd)
cell.colors <- ident.colors[Idents(object = merged_ect.integd)]
names(x = cell.colors) <- colnames(x = merged_ect.integd)
#Prepare to output results as pdf
pdf(file = NULL)
#Parameters
#n = neighborhood size (default = 100)
#scale = velocity scale to use (default = 'log')
#grid.n = number of points to plot for the grid
cc_mds <- 
  show.velocity.on.embedding.cor(emb = Embeddings(object = merged_ect.integd, 
                                                  reduction = "umap"), 
                                 vel = Tool(object = merged_ect.integd, 
                                           slot = "RunVelocity"),
                                 n = 200, 
                                 scale = "sqrt", 
                                 cell.colors = ac(x = cell.colors), 
                                 cex = 0.8,
                                 arrow.scale = 5, 
                                 show.grid.flow = TRUE, 
                                 min.grid.cell.mass = 0.5, 
                                 grid.n = 40,
                                 arrow.lwd = 1, 
                                 do.par = FALSE,
                                 cell.border.alpha = 0.1,
                                 n.cores = 1)
pdf("/publication/supp18_EF_endoderm_500_0221.pdf")
show.velocity.on.embedding.cor(emb = Embeddings(object = merged_ect.integd, 
                                                reduction = "umap"), 
                               vel = Tool(object = merged_ect.integd, 
                                          slot = "RunVelocity"),
                               cc = cc_mds$cc,
                               n = 200, 
                               scale = "sqrt", 
                               cell.colors = ac(x = cell.colors, 
                                                alpha = 0.8), 
                               cex = 0.8,
                               arrow.scale = 5, 
                               show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, 
                               grid.n = 30,
                               arrow.lwd = 1, 
                               do.par = FALSE,
                               cell.border.alpha = 0.1,
                               n.cores = 1)
dev.off()

#Run scVelo
scv <- import("scvelo")
ad <- read_h5ad("/velocity/EF_endoderm_500_0221.h5ad")
ad
ad$obs$idents
#Run scvelo dynamic model
scv$pp$filter_genes(ad) ## filter
scv$pp$moments(ad) ## normalize and compute moments
scv$tl$recover_dynamics(ad) ## model
#Plot (creates pop up window)
scv$tl$velocity(ad, mode='dynamical')
scv$tl$velocity_graph(ad)
ad$obs$idents <- as.factor(ad$obs$idents)
scv$pl$velocity_embedding_stream(ad, 
                                 basis = 'umap', 
                                 color = 'idents',
                                 legend_fontsize = 4,
                                 dpi = 300,
                                 palette = unique(pdata$color),
                                 title = "",
                                 save = "supp18_EF_endoderm_500_0221.png")
ad$obs$stage <- ad$obs$orig.ident
ad$obs$stage <- gsub("E1-1","N4",ad$obs$stage)
ad$obs$stage <- gsub("E1-2","N4",ad$obs$stage)
ad$obs$stage <- gsub("E1-4","N4",ad$obs$stage)
ad$obs$stage <- gsub("E","N4",ad$obs$stage)
ad$obs$stage <- gsub("F1-1","T1",ad$obs$stage)
ad$obs$stage <- gsub("F1-2","T1",ad$obs$stage)
ad$obs$stage <- gsub("F2-1","T1",ad$obs$stage)
ad$obs$stage <- gsub("F2-2","T1",ad$obs$stage)
ad$obs$stage <- gsub("F2-3","T1",ad$obs$stage)
scv$pl$velocity_embedding_stream(ad, 
                                 basis = 'umap', 
                                 color = 'stage',
                                 legend_fontsize = 4,
                                 legend_loc = 'right',
                                 dpi = 300,
                                 title = "",
                                 save = "supp18_EF_endoderm_500_0221_stage.png")
which(ad$obs == "Mesoendoderm", arr.ind=TRUE)
ad$obs$group <- gsub("Mesoendoderm","Endoderm",ad$obs$group)
ad$obs$group <- as.factor(ad$obs$group)
scv$pl$velocity_embedding_stream(ad, 
                                 basis = 'umap', 
                                 color = 'group',
                                 legend_fontsize = 4,
                                 legend_loc = 'right',
                                 dpi = 300,
                                 title = "",
                                 save = "supp18_EF_endoderm_500_0221_group.png")
ad$raw = scv$pp$log1p(ad)
ad$write_h5ad("EF_endoderm_500_0221_analyzed.h5ad")
named_df <- as.matrix(ad$uns$velocity_graph)
rownames(named_df) <- ad$obs$idents
colnames(named_df) <- ad$obs$idents
named_result1 <- rowsum(named_df,rownames(named_df))
named_result2 <- t(rowsum(t(named_result1),colnames(named_result1)))
cell_number <- as.data.frame(table(ad$obs$idents))
normalized_result <- named_result2/cell_number$Freq
normalized_result <- as.data.frame(normalized_result)
normalized_result$names <- cell_number$Var1
result_parsed <- normalized_result[1:24,25:57]
result_parsed$names <- normalized_result$names[1:24]
result_parsed2 <- normalized_result[1:85,86:178]
result_parsed2$names <- normalized_result$names[1:85]
write.csv(result_parsed, row.names = FALSE,
          "/publication/supp18_table_N4_T1_endoderm.csv")
result_parsed <- 
  read.csv("/publication/supp18_table_N4_T1_endoderm.csv", header = T)
library(data.table)
result_parsed_long2 <- melt(setDT(result_parsed), id.vars = c("names"), 
                            variable.name = "B_celltype")
result_parsed_split2 <- split(result_parsed_long2,
                              result_parsed_long2$B_celltype)
listnames <- unique(result_parsed_long2$B_celltype)
df2 = data.frame()
for (i in listnames) {
  table_i <- as.data.frame(table(result_parsed_split2[[i]]$names))
  table_i$To <- result_parsed_split2[[i]]$B_celltype
  table_i$normalized <- result_parsed_split2[[i]]$value
  table_i$sum <- sum(result_parsed_split2[[i]]$value)
  table_i$edge_weight <- table_i$normalized/table_i$sum
  df2 = rbind(df2, table_i)
}
head(df2)
tail(df2)
df2$To <- gsub("\\.",":",df2$To)
write.csv(df, row.names = FALSE, "/publication/supp18_table_N4_T1_endoderm.csv")
df_filtered <- read.csv("/publication/supp18_method1_N4_T1_results.csv",
                        header = T)
df_parsed <- dplyr::filter(df_filtered, Ratio >= 0.2)
df_parsed$edges <- paste0(df_parsed$N4,df_parsed$T1)
df_2_parsed <- dplyr::filter(df2, edge_weight >= 0.2)
df_2_parsed$edges <- paste0(df_2_parsed$Var1,df_2_parsed$To)
write.csv(intersect(df_parsed$edges,df_2_parsed$edges),
          "/publication/supp18_intersect_N4_T1_endoderm_results.csv")
data_new1 <- df[order(df$edge_weight, decreasing = TRUE), ]  # Order data descending
data_new1 <- Reduce(rbind, by(data_new1, data_new1["To"], head, n = 5))
#Top dynamic genes
topgenes <- ad$var["fit_likelihood"]
topgenes_vals <- topgenes[,1]
names(topgenes_vals) <- rownames(topgenes)
topgenes_vals <- sort(topgenes_vals, decreasing=TRUE)
head(topgenes_vals)
write.csv(topgenes_vals, "/velocity/D_all_top_velocity.csv")
scv$pl$scatter(ad, color = 'idents',
               basis = names(topgenes_vals)[1:5], ncols=5, frameon=FALSE,
               save = "/velocity/D_all_velocity_1.png")
import anndata
import pandas as pd
import numpy as np
adata=anndata.read("EF_all_200_0129_analyzed.h5ad")
result = adata.uns['velocity_graph']
B = result.toarray()
i, j = np.triu_indices_from(B, k=1)
v = B[i, j]
ijv = np.concatenate((i, j, v)).reshape(3, -1).T
ijv = ijv[v != 0.0]
df_ijv = pd.DataFrame(ijv)
print(df_ijv)
a_file = open("/velocity/DE_ectoderm_velocity.txt", "w")
np.savetxt(a_file, df_ijv)
a_file.close()
adata.obs[['class','idents']].to_csv("/velocity/DE_ectoderm_identity.csv")
vel <- read.table("/velocity/DE_ectoderm_velocity.txt", header = F)
head(vel)
newdata <- vel[min(which(vel$V1 == 46289)) : max(which(vel$V1 == 46354)) ,]

#k5 neighbor analysis
cluster <- read.csv("/publication/green_sea_urchin/GSU20_GSU24_k5_new_2_5.csv",
           header = T)
cluster_m <- melt(cluster, id = c("GSU24"))
cluster_s <- split(cluster_m, cluster_m$GSU24)
listnames <- as.character(unique(cluster_m$GSU24))
df = data.frame()
for (i in listnames) {
  table_i <- as.data.frame(table(cluster_s[[i]]$value))
  table_i$sample <- i
  table_i$sum <- nrow(cluster_s[[i]])
  table_i$ratio <- table_i$Freq/table_i$sum
  df = rbind(df, table_i)
}
head(df)
tail(df)
write.csv(df, row.names=FALSE,
          "/publication/green_sea_urchin/GSU20_GSU24_k5_new_2_5_analyzed.csv")

#Draw heatmap
colors <- colorRampPalette(c("red","blue"))(50)
gene_fold <- read.csv("/18_samples/pit_tfs_filtered.csv", header = T)
head(gene_fold)
dim(gene_fold)
g_x <- gene_fold[,6:13]
head(g_x)
g_x <- log2(g_x+1)
head(g_x)
rownames(g_x) <- gene_fold[,4]
x <- as.matrix(g_x)
head(x)
pheatmap(x, scale = "none", show_rownames = T, cluster_cols = F,
         cellwidth = 15, cellheight = 3, cutree_col = 1, fontsize = 3,
         legend = FALSE, filename = "pit_markers_filtered.tiff")