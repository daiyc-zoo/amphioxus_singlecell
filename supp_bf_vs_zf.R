### orthofinder:
### Emms, D.M., Kelly, S. 
### OrthoFinder: phylogenetic orthology inference for comparative genomics. 
### Genome Biol 20, 238 (2019).
### https://doi.org/10.1186/s13059-019-1832-y

### code for calculating gene expression homology across species:
### Qiu, C., Cao, J., Martin, B.K. et al. 
### Systematic reconstruction of cellular trajectories across mouse embryogenesis. 
### Nat Genet 54, 328–341 (2022). 
### https://doi.org/10.1038/s41588-022-01018-x
### adopted from https://github.com/ChengxiangQiu/tome_code

#clear objects in workspace
rm(list=ls())
options(stringsAsFactors = F)

#required packages
library("Seurat")
library("dplyr")
library("ggplot2")
library(cowplot)
library("nnls")
library("tidyr")
library(data.table)

#read stored amphioxus rds data
E_all <- readRDS("/publication/analyzed/N4.rds")
Idents(E_all) <- "celltype"
ave <- AverageExpression(object = E_all, return.seurat = TRUE)
write.csv(GetAssayData(object = ave, slot = 'count'),
          "/publication/supp_web_N4.csv")
#read downloaded zebrafish rds data
Z_18 <- readRDS("/singleron/zebrafish/seurat_object_hpf18.rds")
Z_18 <- NormalizeData(Z_18, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
Idents(Z_18) <- "cell_state"
ave <- AverageExpression(object = Z_18, return.seurat = TRUE)
write.csv(GetAssayData(object = ave, slot = 'count'),
          "/publication/supp_zebrafish_18.csv")
#read orthofinder results
ortho <- read.table("/Results_Mar11_1/Orthogroups/Orthogroups.tsv",
                    sep = '\t', header = TRUE)
new_ortho <- cbind(ortho[,1:2],ortho[,4])
colnames(new_ortho) <- c("orthogroup","amphioxus","zebrafish")
new_ortho_clean1 <- new_ortho[!(is.na(new_ortho$amphioxus) | new_ortho$amphioxus==""), ]
new_ortho_clean2 <- new_ortho_clean1[!(is.na(new_ortho_clean1$zebrafish) | new_ortho_clean1$zebrafish==""), ]
new_ortho_clean3 <- setDT(new_ortho_clean2)[, lapply(.SD, function(x) unlist(tstrsplit(x, ",", fixed=TRUE))), by = zebrafish][!is.na(amphioxus)]
new_ortho_clean_final <- setDT(new_ortho_clean3)[, lapply(.SD, function(x) unlist(tstrsplit(x, ",", fixed=TRUE))), by = amphioxus][!is.na(zebrafish)]
new_ortho_clean_final$amphioxus <- gsub("-R.*", "", new_ortho_clean_final$amphioxus)
new_ortho_clean_final$amphioxus <- gsub("_", "-", new_ortho_clean_final$amphioxus)
new_ortho_clean_final$zebrafish <- gsub("\\..*", "", new_ortho_clean_final$zebrafish)
dim(new_ortho_clean_final)
artificial <- paste0("bf_vs_zf_",seq(1:33990))
new_ortho_clean_final <- cbind(new_ortho_clean_final, artificial)
write.csv(new_ortho_clean_final, "/publication/supp_bf_vs_zf_genes.csv")
E_express <- read.csv("/publication/supp_web_N4.csv", header = T)
E_express = E_express[E_express$X %in% new_ortho_clean_final$amphioxus, ]
names(E_express)[1] = "amphioxus"
E_merge <- left_join(E_express,new_ortho_clean_final, by = "amphioxus")
Z_express <- read.csv("/publication/supp_zebrafish_18.csv", header = T)
Z_express= Z_express[Z_express$X %in% new_ortho_clean_final$zebrafish, ]
names(Z_express)[1] = "zebrafish"
Z_merge <- left_join(Z_express,new_ortho_clean_final, by = "zebrafish")
dim(Z_merge)
head(Z_merge)
E_shared <- E_merge[E_merge$artificial %in% Z_merge$artificial, ]
Z_shared <- Z_merge[Z_merge$artificial %in% E_merge$artificial, ]
E_express_new <- E_shared[,-c(1,87:89)]
rownames(E_express_new) <- E_shared[,89]
Z_express_new <- Z_shared[,-c(1,35:37)]
rownames(Z_express_new) <- Z_shared[,37]
E_express_new = E_express_new[ order(row.names(E_express_new)), ]
Z_express_new = Z_express_new[ order(row.names(Z_express_new)), ]
correlation_analysis_nnls_spec_gene <- function(MM1, MM2, fold.change = 1.5, top_gene_num = 300, spec_gene_num = 300) {
  cell_names = colnames(MM1)
  gene_median_expr = apply(MM1, 1, median)
  gene_max_expr = apply(MM1, 1, max)
  correlation_matrix = lapply(cell_names, function(x) {
    gene_other_max_expr = apply(MM1[, colnames(MM1) != x],1,  max)
    target_expr = MM1[, colnames(MM1) == x]
    df_tmp = data.frame(target = target_expr, ratio = (target_expr + 1) / (gene_median_expr + 1), gene_id = row.names(MM1))
    gene_markers = (df_tmp %>% filter(ratio > fold.change) %>%  arrange(desc(ratio)) %>% head(top_gene_num))$gene_id
    df_tmp = data.frame(target = target_expr, ratio = (target_expr + 1) / (gene_other_max_expr + 1), gene_id = row.names(MM1))
    top_markers = (df_tmp %>% arrange(desc(ratio)) %>% head(spec_gene_num))$gene_id
    selected_row = (row.names(MM1) %in% c(as.character(top_markers), as.character(gene_markers)))
    MM1_filtered = MM1[selected_row, ]
    MM2_filtered = MM2[selected_row, ]
    nnls_result = nnls::nnls(MM2_filtered, MM1_filtered[, x])
    df_coef = data.frame("cell_name_MM2" = colnames(MM2), "beta" = coef(nnls_result), "cell_name_MM1" = x)
    return(df_coef)
  })
  result = do.call(rbind, correlation_matrix)
  result = result %>% spread(cell_name_MM2, beta)
  result_2 = result %>% select(-cell_name_MM1)
  rownames(result_2) = result$cell_name_MM1
  return(result_2)
}
correlation_analysis_bidirection <- function(MM1, MM2, fold.change = 1.5, top_gene_num = 300, spec_gene_num = 300) {
  com_1 = correlation_analysis_nnls_spec_gene(MM1, MM2, fold.change = fold.change, top_gene_num = top_gene_num, spec_gene_num = spec_gene_num)
  result = as.data.frame(com_1)
  result$source = rownames(result)
  result_1 = result %>% gather(key = target, value = beta_1, 1:ncol(MM2))
  com_2 = correlation_analysis_nnls_spec_gene(MM2, MM1, fold.change = fold.change, top_gene_num = top_gene_num, spec_gene_num = spec_gene_num)
  result = as.data.frame(com_2)
  result$target = rownames(result)
  result_2 = result %>% gather(key = source, value = beta_2, 1:ncol(MM1))
  result = inner_join(result_1, result_2)
  return(result)
}
conn <- correlation_analysis_bidirection(as.matrix(E_express_new), 
                                         as.matrix(Z_express_new), 
                                         fold.change = 1.5, 
                                         top_gene_num = 1200, 
                                         spec_gene_num = 1200)
conn <- correlation_analysis_bidirection(as.matrix(Z_express_new), 
                                         as.matrix(E_express_new), 
                                         fold.change = 1.5, 
                                         top_gene_num = 1200, 
                                         spec_gene_num = 1200)
library("data.table")
conn$beta <- 2*(conn$beta_1+0.001)*(conn$beta_2+0.001)
dat <- conn[,c("source","target","beta")]
dat <- dcast(dat, source~target)
rownames(dat) <- dat[,1]; dat <- dat[,-1]
write.csv(conn[conn$beta != 2e-6,], "/publication/supp_zf_vs_bf_beta.csv")
write.csv(dat, "/publication/supp_zf_vs_bf_dat.csv")