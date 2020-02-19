#'
#' load library 
#'
library(BiocParallel)
multi_param <- MulticoreParam(worker=2)
register(multi_param, default=TRUE)

library(DESeq2)
library(xlsx)
library(ggplot2)
library(tidyverse)
library(pheatmap)
pkg_dir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsy.all"
fig_dir <- file.path(pkg_dir, "manuscript", "figures")
table_dir <- file.path(pkg_dir, "manuscript", "tables")

#'
#' load common data: time_1, time_1, sanitized.dds, sanitized.rlg
#'                   year2.rlg, year2.dds
#'
source(file.path(pkg_dir, "scripts", "tools.R"))
source(file.path(pkg_dir, "scripts", "manuscript_tools.R"))

load(file.path(pkg_dir, "public_data", "mri_pathology.rda"))
load(file.path(pkg_dir, "public_data", "sanitized.dds.rda"))
load(file.path(pkg_dir, "public_data", "sanitized.rlg.rda"))
load(file.path(pkg_dir, "public_data", "cc_sample.rda"))
load(file.path(pkg_dir, "public_data", "mri_pathology.rda"))

mk <- get(load(file.path(pkg_dir, "data", "FSHD_markers.rda")))$gene_name
mk <- mk[mk %in% rowData(sanitized.dds)$gene_name]
# sanitized.rlg <- rlg(santized.dds, blind=TRUE)

#'
#' (1) PCA of all, colored by pheno_type
#' 
library(ggrepel)
data <- plotPCA(sanitized.rlg, 
                intgroup=c("pheno_type"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
gg <- ggplot(data, aes(PC1, PC2, color=pheno_type)) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #geom_text(label=rownames(data), size=4,
  #          position=position_jitter(width=1, height=1),
  #          hjust=0.5, vjust=1.25, show.legend = F) +
  labs(title="Initial and follow-up visit RNA-seq sample space by PCA") +
  geom_text_repel(aes(label=rownames(data)), size=2.5, show.legend=FALSE) +
  theme_bw() +
  scale_color_brewer(palette="Set2") +
  theme(legend.position="right", 
        plot.title = element_text(hjust = 0.5, size=10))
pdf(file.path(fig_dir, "comb_PCA.pdf"), width=8, height=6)
plot(gg)
dev.off()        

#'
#' (2) Define CC+/DUX4+ samples and run DESeq2 to get markers 
#' 

#'
#' (2a) DESeq2 and GO (.do_goseq) for CC+/DUX4+ vs Controls
#'

#' cc+/dux4+ only, both year i and ii, no cc++
ccp <- cc_sample %>% filter(DUX4_group == "DUX4+") %>%
  filter(cc_intensity != "cc++")
keep <- colnames(sanitized.dds) %in% ccp$sample_name |
  sanitized.dds$pheno_type == "Control"
dux4_dds <- sanitized.dds[, keep]
design(dux4_dds) <- ~ pheno_type
dux4_dds <- DESeq(dux4_dds)
dux4_res <- results(dux4_dds, alpha=0.05, lfcThreshold=2)
summary(dux4_res) # 683 up-regulated

#' tidy up the up-regulated
df_up <- as.data.frame(dux4_res) %>%
  rownames_to_column(var="gencode_id") %>%
  dplyr::filter(padj < 0.05 & log2FoldChange > 0) %>%
  mutate(gene_name = rowData(dux4_dds[gencode_id])$gene_name) %>%
  mutate(gene_id =
           sapply(strsplit(gencode_id, ".", fixed=TRUE), "[[", 1))

df_down <- as.data.frame(dux4_res) %>%
  rownames_to_column(var="gencode_id") %>%
  dplyr::filter(padj < 0.05 & log2FoldChange < 0) %>%
  mutate(gene_name = rowData(dux4_dds[gencode_id])$gene_name) %>%
  mutate(gene_id =
           sapply(strsplit(gencode_id, ".", fixed=TRUE), "[[", 1))

write.csv(df_up, file=file.path(table_dir, "DESeq_CC+DUX4+_vs_Controls_up.csv"))

#' GO analysis
universal <- sapply(strsplit(rownames(dux4_res), ".", fixed=TRUE),
                    "[[", 1)
enriched <- .do_goseq(universal = universal, 
                      selected_genes = df_up$gene_id,
                      return.DEInCat=FALSE, dds=sanitized.dds)
#enriched$DEInCat <- unlist(enriched$DEInCat)                      
write.csv(enriched, file=file.path(table_dir, "GOSeq_CC+DUX4+_vs_Controls.csv"))    

#' select markers that represent each of the enriched biological processes 
#' for later usage!
robust_markers <- list(dux4 = c("ZSCAN4", "KHDC1L", "LEUTX", 
                          "PRAMEF2", "TRIM43"),
                 cellular_matrix=c("PLA2G2A", "COL19A1", "COL1A1",
                                   "COMP", "SFRP2"),
                 cell_cycle=c("CCNA1", "CDKN1A", "CDKN2A", "TNMD"),
                 inflam=c("CCL18", "CCL13", "CCL19", "C6", "C7"), # CCL18 is good for seperate control & FSHD
                 immunoglobulin = c("IGHG4", "IGHGP", "IGHG1"), 
                 muscle = c("PAX7", "PAX3", "MYOG", "MYOM2")) # remove IGLC2

## Note:  IGKV3-11 and IGHA1 are only highly expressed in 01-0029b and 01-0027b

#'
#' (3) Use df_up genes, k-means cluster on PCA of features, visualization
#' 

#' (3a) define markers: df_up + down-regulated in CC+/DUX4- sampels
muscle <- c(
  "MYF5", "MYF6", "MYH1", # newly added
  "MYH13", "MYH2", "MYH4", "MYH6", "MYH7", # newly added
  "MYHAS", "MYL1", "MYL2", #newly added
  "MYO18B", "MYOD1", "MYOG", "MYOM2", "MYOM3",
  "MYOT", "MYOZ1", "MYOZ2", "MYOZ3", "PAX7", "PAX3")
de_markers <- c(df_up$gencode_id, unname(get_ensembl(muscle, sanitized.rlg)))


#' anchar_sample: 01-0041 (control), 32-0002, "01-0029b", "32-0002b1", and "32-0012b"
anchor_sample <- 
  data.frame(label = c("A", "B", "C", "D", "E"),
             sample_name = c("01-0041", "01-0022-1", "01-0029b","32-0002b1", "32-0012b"),
             stringsAsFactors=FALSE)

#' document of feature space and anchor sample for kmean-clustering
tmp <- data.frame(gencode_id = c(unname(get_ensembl(muscle, sanitized.rlg)), df_up$gencode_id), 
                  stringsAsFactors=FALSE) %>%
  mutate(gene_name = rowData(sanitized.dds[gencode_id])$gene_name) %>%
  add_column(feature_type = c(rep("muscle development associated", length(muscle)), 
                              rep("up-regulated genes in CC+ sampels vs control", nrow(df_up))))
file_name <- file.path(table_dir, "suppl_table_4_feature_space_for_kmean_clustering.xlsx")
write.xlsx(tmp, file=file_name, sheetName="Reduced Feature Space", row.names=FALSE)
write.xlsx(anchor_sample, file=file_name, sheetName="Anchor Sample for kmeans", append=TRUE)

#' (3b) k-means (k=5) on PCA of de_markers
k <- 5
data <- assay(sanitized.rlg[de_markers, ])
pca <- prcomp(t(data))
#cl <- kmeans(pca$x, centers=k, iter.max=15, nstart=5)
cl <- kmeans(pca$x, centers=pca$x[anchor_sample$sample_name, ], iter.max=15)

#' (3c) rename cl$cluster
.rename_cluster_5 <- function(sample_name) {
   anchor_samples <- c(A = "01-0041", B = "32-0002", C = "01-0029b",
                       D = "32-0002b1", E = "32-0012b")
   i <- which(anchor_sample$sample_name %in% as.character(sample_name))
   names(anchor_samples)[i]
}

rename_cluster <- function(cl, k) {
  df <- data.frame(sample_name=names(cl$cluster),
                   cluster=cl$cluster)
  tmp <- df %>% group_by(cluster) %>%
    summarise(rename_cluster=.rename_cluster_5(sample_name)) 
  df <- df %>% mutate(renamed=tmp$rename_cluster[cluster])
}

.face_off_cluster <- function(old_name) {
  new_name <- c(A="Mild", B="Moderate", C="IG-High", D="High", E="Muscle-Low")
  factor(new_name[old_name], levels=new_name)
}

cluster_color <- c(A_Cntr="#ff7f00", A="#a65628", B="#f781bf", C="#984ea3", D="#e41a1c", E="#377eb8")
cluster_df <- rename_cluster(cl, k) %>%
  mutate(RNA_cluster = ifelse(sanitized.rlg$pheno_type == "Control",
                              paste0(renamed, "_Cntr"), renamed)) %>%
  mutate(RNA_cluster = factor(RNA_cluster, levels=c("A_Cntr", "A", "B", "C", "D", "E"))) %>%
  mutate(color = cluster_color[as.character(RNA_cluster)])  %>%
  mutate(new_cluster_name = .face_off_cluster(renamed)) 

save(cluster_df, file=file.path(pkg_dir, "public_data", "cluster_df.rda"))
save(cluster_df, file=file.path(pkg_dir, "data", "cluster_df.rda"))
sanitized.rlg$RNA_cluster_5 <- sanitized.dds$RNA_cluster_5 <- cluster_df$RNA_cluster
#save(sanitized.dds, file=file.path(pkg_dir, "data", "updated_sanitized.dds.rda"))
#save(sanitized.rlg, file=file.path(pkg_dir, "data", "updated_sanitized.rlg.rda"))
 
 
#df <- data.frame(sample_name=names(cl$cluster), cluster=cl$cluster)
# sanitized.rlg$RNA_cluster_5 <- factor(df$cluster)

#' (3d) visualization: PCA showing k-means clustering
library(ggrepel)
data <- plotPCA(sanitized.rlg, 
                intgroup=c("RNA_cluster_5", "pheno_type"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
gg <- ggplot(data, aes(PC1, PC2, color=RNA_cluster_5)) +
  geom_point(size=0.8) +
  #xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  #ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  labs(title="Principal components on initial and follow-up visits") +
  geom_text_repel(aes(label=rownames(data)), size=2.2, show.legend=FALSE) +
  theme_bw() +
  scale_color_manual(name="Classes", values=cluster_color) +
  #scale_color_brewer(palette="Set2") +
  theme(legend.position="right", 
        plot.title = element_text(hjust = 0.5, size=10)) 
    
#  scale_color_discrete(name = "Cluster")

pdf(file.path(fig_dir, "comb_kmeans_cluster_pca.pdf"), width=6, height=4.5)
plot(gg)
dev.off()     

#'
#' (3e) viz: heatmap of all de_markers (supplimental) and robust_markers (figure) with pathology data
#'
robust_markers <- list(dux4 = c("ZSCAN4", "KHDC1L", "LEUTX", 
                          "PRAMEF2", "TRIM43"),
                 cellular_matrix=c("PLA2G2A", "COL19A1", "COL1A1",
                                   "COMP", "SFRP2"),
                 cell_cycle=c("CCNA1", "CDKN1A", "CDKN2A"),
                 inflam=c("CCL18", "CCL13", "C1QA", "C6", "C7"), # CCL18 is good for seperate control & FSHD
                 immunoglobulin = c("IGHG4", "IGHGP", "IGHG1"), 
                 muscle = c("PAX7", "PAX3", "MYOG", "MYOD1")) 
library(pheatmap)

#' data
marker_id <- unlist(lapply(robust_markers, get_ensembl, sanitized.rlg))
data_lst <- lapply(levels(sanitized.rlg$RNA_cluster_5), function(x, rlg) {
  s <- colnames(rlg)[rlg$RNA_cluster_5 == x]
  data <- assay(rlg[marker_id, s]) 
  # re-arrange by the colmeans (look heatmap)
  Colv <- colMeans(data, na.rm = TRUE)
  data <- data[, order(Colv)]
}, rlg=sanitized.rlg)
data <- do.call(cbind, data_lst)
rownames(data) <- rowData(sanitized.rlg[rownames(data)])$gene_name
zscore_data <- (data - rowMeans(data)) / rowSds(data)
gaps_row <- cumsum(lapply(robust_markers, length))

#' annotation bar
classes <- as.character(sanitized.rlg[, colnames(data)]$RNA_cluster_5)
pheno_type <- sanitized.dds[, colnames(data)]$pheno_type                         
annotation_col <- 
  data.frame(sample_name = colnames(data),
             classes     = classes,
             k_means     = sanitized.dds[, colnames(data)]$RNA_cluster_5,
             FSHD = sanitized.dds[, colnames(data)]$pheno_type) 

new_time <- bind_rows(mri_pathology$time_1, mri_pathology$time_2) %>% 
  dplyr::filter(sample_name!="32-0008b") %>%
  mutate(sample_name=as.character(sample_name))

new_time$sample_name[new_time$sample_name=="32-0002b"] <- "32-0002b1"
gaps_col <- cumsum(count(annotation_col, k_means)[, "n"])$n

annotation_col <- annotation_col %>%
  left_join(new_time, by="sample_name") %>%
  rename(Pathology=Pathology.Score, T1=fat_fraction) %>%
  dplyr::select(sample_name, classes, FSHD, T1, STIR_rating, Inflammation, Pathology) %>%
  column_to_rownames(var="sample_name") 

col <- c("#f7f7f7", "#238b45")
ann_cor <- list(
  STIR_rating= col,  T1 = col, FSHD=c(Control="white", FSHD="seagreen"),
  Inflammation = col, Fibrosis = col,
  Var.in.Fiber.Size = col, Nec.Reg.Inflamm = col,
  Nucleation = col, Pathology = col,
  classes = cluster_color)

#' heatmap of robust markers
pheatmap(zscore_data, annotation_col=annotation_col,
         fontsize_row=7, fontsize_col=6,
         scale="none", 
         gaps_col=gaps_col,
         gaps_row=gaps_row,
         annotation_colors=ann_cor,
         cellheight=10, annotation_legend=FALSE,
         cluster_rows=FALSE, cluster_cols=FALSE,
         file=file.path(fig_dir, "comb_kmeans_cluster_heatmap.pdf"))

#' 
#' supplemental figure: heatmap of all markers
#'
data_lst <- lapply(levels(sanitized.rlg$RNA_cluster_5), function(x, rlg) {
  s <- colnames(rlg)[rlg$RNA_cluster_5 == x]
  data <- assay(rlg[de_markers, s]) 
  # re-arrange by the colmeans (look heatmap)
  Colv <- colMeans(data, na.rm = TRUE)
  data <- data[, order(Colv)]
}, rlg=sanitized.rlg)
data <- do.call(cbind, data_lst)
rownames(data) <- rowData(sanitized.rlg[rownames(data)])$gene_name
zscore_data <- (data - rowMeans(data)) / rowSds(data)
annotation_col.2 <- annotation_col[colnames(zscore_data), ]
pheatmap(zscore_data, annotation_col=annotation_col.2,
         fontsize_row=7, fontsize_col=6,
         show_rownames=FALSE,
         scale="none", 
         gaps_col=gaps_col,
         annotation_colors=ann_cor,
         annotation_legend=FALSE,
         cluster_rows=FALSE, cluster_cols=FALSE,
         file=file.path(fig_dir, "comb_kmeans_cluster_heatmap.allDE.pdf"))

#'
#' (4) DESeq for each class of FSHDs
#'
#load(file.path(pkg_dir, "data", "updated_sanitized.dds.rda"))
#load(file.path(pkg_dir, "data", "updated_sanitized.rlg.rda"))

.tidy_results <- function(res, dds, padj=0.05) {
  de <- as.data.frame(res) %>%
    rownames_to_column(var="gencode_id") %>%
    filter(padj < 0.05) %>%
    mutate(gene_name = rowData(dds[gencode_id])$gene_name) %>%
    mutate(gene_id =
             sapply(strsplit(gencode_id, ".", fixed=TRUE), "[[", 1))
}

#'
#' (4a) FSHD class A (mild) vs controls
#' log2FoldChanges > 1 with adj-pvalue < 0.05
#'
new.dds <- sanitized.dds
design(new.dds) <- ~ RNA_cluster_5
new.dds <- DESeq(new.dds)
classA.res <- results(new.dds, alpha = 0.05, 
                      name = "RNA_cluster_5_A_vs_A_Cntr")
de <- .tidy_results(classA.res, dds=new.dds, padj=0.05)                      
de_up <- de %>% dplyr::filter(log2FoldChange > 1) %>%
  mutate(DUX4_marker=gene_name %in% mk)
de_down <- de %>% dplyr::filter(log2FoldChange < 0) 

universal <- sapply(strsplit(rownames(classA.res), ".", fixed=TRUE),
                    "[[", 1)
enriched <- .do_goseq(universal = universal, 
                      selected_genes = de_up$gene_id,
                      return.DEInCat=TRUE, dds=sanitized.dds)

file_name <- file.path(table_dir, "MildFSHD_vs_Controls.xlsx")
write.xlsx(de_up, sheetName="Up_regulated_Genes", file=file_name,
           append=FALSE)
write.xlsx(enriched, sheetName="GOSeq_upDE_enriched", file=file_name,
           append=TRUE)
write.xlsx(de_down, sheetName="Down_regulated_Genes", file=file_name,
           append=TRUE)           


###
#' Figure 5C: condensed plot
###
df_hot <- data.frame(term = c("DUX4", "immune/inflam", "extracellular", "immunoglobulin", "muscle"),
                     intensity = c(4, 4, 4, 2, 4 ),
                     class = "High")
df_mild <- data.frame(term = c("DUX4", "immune/inflam", "extracellular", "immunoglobulin", "muscle"),
                     intensity = c(0.5, 0.5, 1, 0, 4 ),
                     class = "Mild")
df_moderate <- data.frame(term = c("DUX4", "immune/inflam", "extracellular", "immunoglobulin", "muscle"),
                     intensity = c(2, 1, 2, 1, 4 ),
                     class = "Moderate")
df_img <- data.frame(term = c("DUX4", "immune/inflam", "extracellular", "immunoglobulin", "muscle"),
                     intensity = c(2, 3, 3, 4, 3 ),
                     class = "Ig-High")
df_muscle <- data.frame(term = c("DUX4", "immune/inflam", "extracellular", "immunoglobulin", "muscle"),
                     intensity = c(0, 3, 3, 2, 0 ),
                     class = "Muscle-Low")
df_lst <- list(df_hot, df_mild, df_moderate, df_img, df_muscle)
df <- do.call(rbind, df_lst)
df$class <- factor(df$class, levels=c("Mild", "Moderate", "Ig-High", "High", "Muscle-Low"))
library(ggthemes)
gg <- ggplot(df, aes(y=intensity, x=term, group=class)) +
  geom_bar(stat="identity", width=0.4) +
  coord_flip() +
  #theme_wsj() +
  theme_economist() +
  #theme_minimal() +
  facet_wrap(~class, nrow = 1, strip.position = "bottom") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = 9))
pdf(file.path(fig_dir, "comb_cluster_intensity.pdf"), width=5, height=1.5)
gg
dev.off()  
                                                                                                                       
#'
#' MRI and marker scores (adjust_time_1 and time_2?)
#'

load(file.path(pkg_dir, "public_data", "cluster_df.rda"))
load(file.path(pkg_dir, "public_data", "mri_pathology.rda"))
load(file.path(pkg_dir, "public_data", "sanitized.rlg.rda"))
santized.rlg$RNA_cluster_5 <- cluster_df$RNA_cluster
#' fix time_var sample_name: 32-0002b -> 32-0002b1
histo <- bind_rows(mri_pathology$time_1, mri_pathology$time_2) %>%
  mutate(sample_name = as.character(sample_name)) %>%
  mutate(sample_name = ifelse(sample_name == "32-0002b", "32-0002b1", sample_name)) %>%
  left_join(cluster_df, by="sample_name")

#' sanity check: all cluster$sample_name are in time_var$sample_name
all(cluster_df$sample_name %in% histo$sample_name)

#' join cluster_df to histology/mri scores
mri_stats <- histo %>% 
  group_by(new_cluster_name)  %>%
  summarize(pathology_avg = mean(Pathology.Score, na.rm=TRUE),
            pathology_sd = sd(Pathology.Score, na.rm=TRUE),
            T1_fraction_avg = mean(fat_fraction, na.rm=TRUE),
            T1_fraction_sd = sd(fat_fraction, na.rm=TRUE),
            STIR_avg = mean(STIR_rating),
            STIR_sd = sd(STIR_rating),
            Inflam_avg = mean(Inflammation, na.rm=TRUE),
            Infalm_se = mean(Inflammation, na.rm=TRUE))
write.csv(mri_stats, file=file.path(table_dir, "cluster_vs_pathology_and_FF.csv"))
save(mri_stats, file=file.path(pkg_dir, "data", "mri_stats_per_cluster.rda"))

#' visualize cluster vs pathology score and fat fraction
gg_path <- ggplot(histo, aes(x=new_cluster_name, y=Pathology.Score, group=new_cluster_name)) +
  geom_boxplot(width=0.7, outlier.size = 1) + 
  theme_bw() +
  labs(y="Pathology score", title="Pathology score vs. RNA-seq category") +
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(angle=90, hjust=1))

gg_inflam <- ggplot(histo, aes(x=new_cluster_name, y=Inflammation, group=new_cluster_name)) +
  geom_boxplot(width=0.7, outlier.size = 1) + 
  theme_bw() +
  labs(y="Inflammation infiltrates", title="Inflammation vs. RNA-seq category") +
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(angle=90, hjust=1))

gg_T1 <- ggplot(histo, aes(x=new_cluster_name, y=fat_fraction, group=new_cluster_name)) +
  geom_boxplot(width=0.7, outlier.size = 1) + 
  #geom_jitter(size=0.7) +
  theme_bw() +
  labs(y="T1 fraction", title="T1 fraction vs. RNA-seq category") +
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(angle=90, hjust=1))

gg_stir <- ggplot(histo, aes(x=new_cluster_name, y=STIR_rating, group=new_cluster_name)) +
  geom_boxplot(width=0.7, outlier.size = 1) + 
  #geom_jitter(size=0.7) +
  theme_bw() +
  labs(y="STIR rating", title="STIR rating vs. RNA-seq category") +
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(angle=90, hjust=1))

library(gridExtra)
pdf(file.path(fig_dir, "comb_cluster_vs_mri_path.pdf"), width=8, height=4)
grid.arrange(gg_path, gg_inflam, gg_T1, nrow=1)
dev.off()

#'
#' try out:
#'
col <- as.character(wes_palette(n=4, name="Cavalcanti1", type="discrete"))
tmp <- mri_stats %>%
  dplyr::select(new_cluster_name, pathology_avg, Inflam_avg, STIR_avg) %>%
  rename(Pathology=pathology_avg, Inflamm=Inflam_avg, STIR=STIR_avg) %>%
  gather(feature, average.score, -new_cluster_name)
gg <- ggplot(tmp, aes(x=new_cluster_name, y=average.score, fill=feature)) +
  geom_bar(stat="identity", width=0.7)  +
  #scale_fill_brewer(palette="Set2") +
  scale_fill_manual(values=col) +
  theme_bw() +
  labs(y="Average score") +
  theme(plot.title = element_text(hjust = 0.5, size=9),
        legend.key.size = unit(0.8, "line"),
        #legend.position="bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(angle=25, hjust=1))
pdf(file.path(fig_dir, "comb_mri_cluster_stats.pdf"), width=3, height=2.5)
plot(gg)
dev.off()

#'
library(xlsx)
file_name <- file.path(table_dir, "Suppl_table_6_Candidate_Biomarkers_and_Enriched_Go.xlsx")
de <- read.xlsx(file=file_name, sheetIndex=1)
ecm_de <- de %>% filter(ecm == TRUE) %>% pull(gencode_id) %>% as.character(.)
dux4_de <- de %>% filter(DUX4_induced == TRUE) %>% pull(gencode_id) %>% as.character(.)
inflamm_de <- de %>% filter(inflamm == TRUE | immune == TRUE) %>% 
  pull(gencode_id) %>% as.character(.)
stress_de <- de %>% filter(stress == TRUE) %>% pull(gencode_id) %>% as.character(.)
ig_de <- de %>% filter(grepl("IGH", gene_name) | grepl("IGK", gene_name)) %>%
  pull(gencode_id) %>% as.character(.)
de_list <- list(Extracellur_Matrix=ecm_de, DUX4=dux4_de, Inflammatory=inflamm_de, Immunoglobulin=ig_de)

#' add muscle
muscle <- c(
  "MYF5", "MYF6", "MYH1", # newly added
  "MYH13", "MYH2", "MYH4", "MYH6", "MYH7", # newly added
  "MYHAS", "MYL1", "MYL2", #newly added
  "MYO18B", "MYOD1", "MYOG", "MYOM2", "MYOM3",
  "MYOT", "MYOZ1", "MYOZ2", "MYOZ3", "PAX7", "PAX3")
de_list$Muscle <- unname(get_ensembl(muscle, sanitized.rlg))

# average expression relative to controls
de_avg_per_class <- lapply(de_list, function(id){
  tmp <- sanitized.rlg[id, ]
  class <- levels(sanitized.rlg$RNA_cluster_5)
  avg_rlg <- sapply(class, function(x){
    colMeans(assay(tmp)[, tmp$RNA_cluster_5 == x])
  })
  avg_per_class <- sapply(avg_rlg, mean)
  avg_per_class - avg_per_class["A_Cntr"]
})

de_avg_per_class <- as.data.frame(do.call(cbind, de_avg_per_class)) %>%
  rownames_to_column(var="class") %>% 
  dplyr::filter(class != "A_Cntr") %>%
  mutate(RNA_cluster = .face_off_cluster(class)) %>%
  dplyr::select(-class) %>%
  gather(feature, average.score, -RNA_cluster) %>%
  mutate(feature=factor(feature))

col <- as.character(wes_palette(n=5, name="Darjeeling1", type="discrete"))
gg <- ggplot(de_avg_per_class, aes(x=RNA_cluster, y=average.score, fill=feature)) +
  geom_bar(stat="identity", width=0.75)  +
  #scale_fill_brewer(palette="Set2") +
  scale_fill_manual(values=col) +
  theme_bw() +
  labs(y="Average expression - control") +
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=9),
        axis.text.x = element_text(angle=25, hjust=1, size=8), 
        legend.key.size = unit(0.75, "line"))
pdf(file.path(fig_dir, "comb_cluster_avg_marker_scores.pdf"), width=3.2, height=2.5)
plot(gg)
dev.off()

