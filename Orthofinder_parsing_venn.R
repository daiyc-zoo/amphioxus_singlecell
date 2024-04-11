aTF <- read.csv("a_all_tf.csv", header = T)
mTF <- read.csv("m_all_tf.csv", header = T)
zTF <- read.csv("z_all_tf.csv", header = T)
tTF <- read.csv("t_all_tf.csv", header = T)
uTF <- read.csv("u_all_tf.csv", header = T)
nTF <- read.csv("n_all_tf.csv", header = T)
ortho <- 
  read.table(file = '/data/home/slr/dyc/singleron/SAMap/orthofinder/proteomes/OrthoFinder/Results_Sep27/Orthogroups/Orthogroups.tsv', 
             sep = '\t', header = TRUE)
###
library(tidyr)
###
amp_only <- ortho[,c(1,2)]
amp_only_clean <- 
  separate_rows(amp_only,Branchiostoma_floridae_clean2,sep=",\\s+")
amp_only_clean$X <- amp_only_clean$Branchiostoma_floridae_clean2
amp_only_clean$Branchiostoma_floridae_clean2 <- NULL
TF_interest <- aTF$X
amp_filter <- amp_only_clean[amp_only_clean$X %in% TF_interest,]
amp_total <- merge(amp_filter, aTF, by="X", all = TRUE)
write.csv(amp_total, "a_all_tf_ortho_annot.csv")
###
mouse_only <- ortho[,c(1,6)]
mouse_only_clean <- 
  separate_rows(mouse_only,Mus_musculus_clean,sep=",\\s+")
mouse_only_clean$X <- mouse_only_clean$Mus_musculus_clean
mouse_only_clean$Mus_musculus_clean <- NULL
TF_interest <- mTF$X
mouse_filter <- mouse_only_clean[mouse_only_clean$X %in% TF_interest,]
mouse_total <- merge(mouse_filter, mTF, by="X", all = TRUE)
write.csv(mouse_total, "m_all_tf_ortho_annot.csv")
###
z_only <- ortho[,c(1,4)]
z_only_clean <- 
  separate_rows(z_only,Danio_rerio_clean,sep=",\\s+")
z_only_clean$X <- z_only_clean$Danio_rerio_clean
z_only_clean$Danio_rerio_clean <- NULL
TF_interest <- zTF$X
z_filter <- z_only_clean[z_only_clean$X %in% TF_interest,]
z_total <- merge(z_filter, zTF, by="X", all = TRUE)
write.csv(z_total, "z_all_tf_ortho_annot.csv")
###
tTF$X <- gsub("KH2012:","KH2012_",tTF$X)
t_only <- ortho[,c(1,3)]
t_only_clean <- 
  separate_rows(t_only,Ciona_intestinalis_KH2012,sep=",\\s+")
t_only_clean$X <- t_only_clean$Ciona_intestinalis_KH2012
t_only_clean$Ciona_intestinalis_KH2012 <- NULL
TF_interest <- tTF$X
t_filter <- t_only_clean[t_only_clean$X %in% TF_interest,]
t_total <- merge(t_filter, tTF, by="X", all = TRUE)
write.csv(t_total, "t_all_tf_ortho_annot.csv")
###
u_only <- ortho[,c(1,5)]
u_only_clean <- 
  separate_rows(u_only,L_var_proteins_clean3,sep=",\\s+")
u_only_clean$X <- u_only_clean$L_var_proteins_clean3
u_only_clean$L_var_proteins_clean3 <- NULL
TF_interest <- uTF$X
u_filter <- u_only_clean[u_only_clean$X %in% TF_interest,]
u_total <- merge(u_filter, uTF, by="X", all = TRUE)
write.csv(u_total, "u_all_tf_ortho_annot.csv")
###
n_only <- ortho[,c(1,7)]
n_only_clean <- 
  separate_rows(n_only,Nematostella_proteins,sep=",\\s+")
n_only_clean$X <- n_only_clean$Nematostella_proteins
n_only_clean$Nematostella_proteins <- NULL
TF_interest <- nTF$X
n_filter <- n_only_clean[n_only_clean$X %in% TF_interest,]
n_total <- merge(n_filter, nTF, by="X", all = TRUE)
write.csv(n_total, "n_all_tf_ortho_annot.csv")
###
library(venn)
x <- list(Amphioxus = unique(na.omit(amp_total$Orthogroup)), 
          Mouse = unique(na.omit(mouse_total$Orthogroup)), 
          Zebrafish = unique(na.omit(z_total$Orthogroup)),
          Tunicate = unique(na.omit(t_total$Orthogroup)),
          Sea_urchin = unique(na.omit(u_total$Orthogroup)),
          Sea_anemone = unique(na.omit(n_total$Orthogroup)))
jpeg("venn_fig5_orthologues.jpg", width = 1000, height = 1000)
venn(x, ilabels = TRUE, zcolor = "style", ilcs = 3, sncs = 3)
dev.off()
###
x <- list(Amphioxus = unique(na.omit(amp_total$Orthogroup)), 
          Mouse = unique(na.omit(mouse_total$Orthogroup)), 
          Zebrafish = unique(na.omit(z_total$Orthogroup)),
          Tunicate = unique(na.omit(t_total$Orthogroup)),
          Sea_urchin = unique(na.omit(u_total$Orthogroup)))
jpeg("venn_fig5_deuterostomes.jpg", width = 1000, height = 1000)
venn(x, ilabels = TRUE, zcolor = "style", ilcs = 3, sncs = 3)
dev.off()
###
x <- list(Amphioxus = unique(na.omit(amp_total$Orthogroup)), 
          Mouse = unique(na.omit(mouse_total$Orthogroup)), 
          Zebrafish = unique(na.omit(z_total$Orthogroup)),
          Tunicate = unique(na.omit(t_total$Orthogroup)))
jpeg("venn_fig5_chordates.jpg", width = 1000, height = 1000)
venn(x, ilabels = TRUE, zcolor = "style", ilcs = 3, sncs = 3)
dev.off()
###
x <- list(Amphioxus = unique(na.omit(amp_total$Orthogroup)), 
          Mouse = unique(na.omit(mouse_total$Orthogroup)), 
          Zebrafish = unique(na.omit(z_total$Orthogroup)))
jpeg("venn_fig5_mouse_zebrafish_amphioxus.jpg", width = 1000, height = 1000)
venn(x, ilabels = TRUE, zcolor = "style", ilcs = 3, sncs = 3)
dev.off()
###
x <- list(Amphioxus = unique(na.omit(amp_total$Orthogroup)), 
          Tunicate = unique(na.omit(t_total$Orthogroup)),
          Sea_urchin = unique(na.omit(u_total$Orthogroup)),
          Sea_anemone = unique(na.omit(n_total$Orthogroup)))
jpeg("venn_fig5_invertebrate.jpg", width = 1000, height = 1000)
venn(x, ilabels = TRUE, zcolor = "style", ilcs = 3, sncs = 3)
dev.off()
###
x <- list(Amphioxus = unique(na.omit(amp_total$Orthogroup)), 
          Tunicate = unique(na.omit(t_total$Orthogroup)))
jpeg("venn_fig5_protochordate.jpg", width = 1000, height = 1000)
venn(x, ilabels = TRUE, zcolor = "style", ilcs = 3, sncs = 3)
dev.off()
###
Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}
Intersect(x)
