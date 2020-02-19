#' this script includes the frist part of RNA-seq analysis for manuscript.

#' NOTE: for year II samples analysis, the normalization only uses year II samples and some 
#' controls with comparable library size from year I batches. 
#' Control "2578", "2550", "2559" are excluded if year II analysis because the libraray sizes
#' are too small, which might result in larger rlog score and I did not like it.

#' NOTE: 32-0008b has too much fat content and therefore exluded from the analysis
#' NOTE: 32-0010b, 32-0012b and 32-0016b have very little muscle content and perhaps that's
#' why DUX4 is not present. 

#'
#' load library
#'
library(DESeq2)
library(xlsx)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(BiocParallel)
multi_param <- MulticoreParam(worker=2)
register(multi_param, default=TRUE)

pkg_dir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsy.all"
fig_dir <- file.path(pkg_dir, "manuscript", "figures")
table_dir <- file.path(pkg_dir, "manuscript", "tables")
sup_dir <- file.path(pkg_dir, "manuscript", "sup_figs")

#' Make Table 2


#'
#' load common data: time_1, time_1, sanitized.dds, sanitized.rlg
#'                   year2.rlg, year2.dds, year2.rlg, year2.dds
#'
source(file.path(pkg_dir, "scripts", "tools.R"))
source(file.path(pkg_dir, "scripts", "manuscript_tools.R"))

load(file.path(pkg_dir, "public_data", "sanitized.dds.rda"))
load(file.path(pkg_dir, "public_data", "mri_pathology.rda"))
load(file.path(pkg_dir, "public_data", "year2.dds.rda"))
load(file.path(pkg_dir, "public_data", "year2.rlg.rda"))

#' Only include follow-up visit dataset and normalize with controls
#' keep these codes for gitbook
#keep_column_2 <- (sanitized.dds$visit == "II" | sanitized.dds$pheno_type == "Control") 

#' keep controls, exclude 32-0008b and 2578  
#year2.dds <- sanitized.dds[, keep_column_2]
#year2.dds$batch <- factor(year2.dds$batch)
#' re-normalized again
#year2.dds <- estimateSizeFactors(year2.dds)
#year2.dds <- estimateDispersions(year2.dds)
#design(year2.dds) <- ~ pheno_type 
#year2.rlg <- rlog(year2.dds, blind = FALSE) #normalized with respect to year2 library size

#'
#' (1) marker scores:rlog sum of markers; append to year2.dds; output excel sheet; 
#' dux4.score and group (hclust)
#' 
markers_list <- list(dux4 = c("LEUTX", "KHDC1L", "PRAMEF2", "TRIM43"),
                     extracellular_matrix = c("PLA2G2A", "COL19A1", "COMP", "COL1A1"), # remove SFPR2
                     inflamm = c("CCL18", "CCL13" ,"C6", "C7"), # remove CCL19
                     cell_cycle = c("CCNA1", "CDKN1A", "CDKN2A"), # newly added
                     immunoglobulin = c("IGHA1", "IGHG4", "IGHGP")) #delete IGKV3-11
#' note: immuno_globin candidate can also be     
#' c("IGHA1", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGKV3-11", "IGHGP", "IGLC2")                                  
markers_id_list <- lapply(markers_list, get_ensembl, rse=year2.dds)
relative_rlg <- lapply(markers_id_list, .get_sum_relative_local_rlg,
                       dds = year2.dds)  

#' (1a) Append markers score to colData of year2.dds and year2.rlg
df <- as.data.frame(do.call(cbind, relative_rlg)) %>%
  dplyr::rename(dux4.rlogsum=dux4, ecm.rlogsum=extracellular_matrix, 
                inflamm.rlogsum=inflamm, cellcycle.rlogsum=cell_cycle,
                img.rlogsum=immunoglobulin)       

if (!any(names(colData(year2.dds)) %in% "dux4.rlogsum")) {
  colData(year2.dds) <- append(colData(year2.dds), 
                               as(df[colnames(year2.dds), ], "DataFrame"))
  colData(year2.rlg) <- colData(year2.dds)
}            


#' (1b) Excel sheet: markers score, count, normalized count, rlog, TPM, RPKM
df <- df %>% 
  add_column(pheno_type=year2.dds[, rownames(df)]$pheno_type) %>%
  rownames_to_column(var="sample_name") %>%
  add_column(non.dux4.avg=rowMeans(df[, c("ecm.rlogsum", "inflamm.rlogsum", 
                                          "cellcycle.rlogsum", "img.rlogsum")])) %>%
  arrange(dux4.rlogsum) #arrange(non.dux4.avg)

file_name <- file.path(table_dir, "Year2_markers_score.xlsx")
write.xlsx(df, sheetName="Markers relative log2 sum",
           row.names=FALSE,
           file=file_name, append=FALSE)
#' counts/normalized counts/TPM/RPKM           
sub_dds <- year2.dds[unlist(markers_id_list)]          
.write_expression_to_xlsx(sub_dds=sub_dds, file_name=file_name,
                          append=TRUE) 
#' local rlog
local_rlg <- lapply(markers_id_list, .get_local_rlg, year2.dds)
local_rlg <- .insert_gene_name(do.call(rbind, local_rlg), year2.dds)
write.xlsx(local_rlg, sheetName="Local rlog", # Blind=TRUE
           file=file_name, append=TRUE)        

#' (1c) duxs.score and dux4.group; append to year2.dds
dux4_local_rlog <- .get_local_rlg(markers_id_list$dux4, year2.dds)
# sanity check
all(rownames(dux4_local_rlog) == markers_id_list$dux4)
dux4_tree <- .get_DUX4_clust_tree(dux4_rlog=dux4_local_rlog)
year2.dds$dux4.group <- year2.rlg$dux4.group <- dux4_tree[colnames(year2.dds)] 

#' (1d) append dux4.group and counts to markers_score.xlsx
file_name <- file.path(table_dir, "Year2_markers_score.xlsx")
.makeFourBiomarkerGroupTable(dds=year2.dds, 
                             relative_logsum = relative_rlg$dux4,
                             dux4_group = dux4_tree,
                             markers_id = markers_id_list$dux4, 
                             file_name = file_name)

#' save some objects for part 2 analysis: year2.dds and year2.rlg
year2_dux4_score <- colData(year2.dds)[, c("sample_name", "dux4.rlogsum", "dux4.group")]
save(year2_dux4_score, file=file.path(pkg_dir, "public_data", "year2_dux4_score.rda"))
save(year2.dds, file=file.path(pkg_dir, "data", "year2.dds.rda"))


#' (1e) heatmap of 16 markers
library(pheatmap)
library(reshape2)
tmp <- markers_id_list
alt_dux4 <- get_ensembl(c("ZSCAN4", "MBD3L2", "PRAMEF25"), year2.dds)
tmp$dux4 <- alt_dux4
#tmp$dux4 <- c(tmp$dux4, alt_dux4)
selected <- melt(tmp)
#data <- assay(year2.rlg[as.character(selected$value)])
#rownames(data) <- rowData(year2.rlg[rownames(data)])$gene_name
#' revise data
data <- as.matrix(.get_local_rlg(as.character(selected$value), year2.dds))
all(rownames(data)==as.character(selected$value))
data <- data[as.character(selected$value), ]
rownames(data) <- rowData(year2.dds[rownames(data)])$gene_name
zscore_data <- (data - rowMeans(data)) / rowSds(data)

#' end of revise data


annotation_col <- data.frame(pheno_type=year2.dds$pheno_type, 
                             dux4_group=factor(year2.dds$dux4.group))
rownames(annotation_col) <- colnames(year2.dds)
annotation_row <- data.frame(marker_type=selected$L1)
rownames(annotation_row) <- rownames(data)
row_gaps <- cumsum(sapply(tmp, length))
row_gaps <- cumsum(sapply(tmp, length))

pheatmap(zscore_data, annotation_col=annotation_col, cellheight=8,
         fontsize_row = 7, fontsize = 7,
         cluster_rows = FALSE, gaps_row = row_gaps,
         annotation_row=annotation_row, scale="none", silent=TRUE,
         file=file.path(fig_dir, "Year2_16Biomarkers_heatmap.pdf"))


#'
#' (2) DUX4 candidate biomarker elevated in FSHD samples
#'

#' (2a) scatter plot of four DUX4 markers
file_name <- file.path(fig_dir, "Year2_Scatter_FourBiomarkers.pdf")
source(file.path(pkg_dir, "scripts", "manuscript_tools.R"))
.makeFourBiomarkerScatterPlot(dds=year2.dds, 
                              group=year2.dds$dux4.group,
                              markers_id=markers_id_list$dux4,
                              call_threshold=2,
                              file_name=file_name,
                              plot.tile_top=FALSE)

#' (2b) supplemental heatmap
.makeFourBiomarkerHeatmap(rlg=year2.rlg,
                          group=year2.rlg$dux4.group,
                          markers_id=markers_id_list$dux4,
                          file_name=file.path(fig_dir, 
                                              "Year2_Heatmap_FourBiormarkers.pdf"))

#' (2c) supplementary heatmap for 53 markers
mk <- get(load(file.path(pkg_dir, "data", "FSHD_markers.rda")))

#'
#' (3) Other biomarkers elevated in FSHD samples by loading vectors of PCA
#'

#' (3a) PCA 
library(ggrepel)
data <- plotPCA(year2.rlg, 
                intgroup=c("pheno_type"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
gg <- ggplot(data, aes(PC1, PC2, color=pheno_type)) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  labs(title="PCs on the follow-up RNA-seq muscle biopsies") +
  geom_text_repel(aes(label=rownames(data)), size=3, show.legend=FALSE) +
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(#legend.position="bottom", 
        legend.title = element_blank(),
        legend.position = c(0.1, 0.1),
        legend.box.background = element_rect(colour = "black"),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size=10))
pdf(file.path(fig_dir, "Year2_PCA.pdf"), width=5, height=4)
plot(gg)
dev.off()  

#' (3b) boxplot of markers (dux4, ecm, cc, img) using TPM
markers_list$candidates <- c("CDKN1A", "CCL13", "COL19A1", "COMP")
markers_exp_list <- lapply(markers_list, .getExpFromDDS, dds=year2.dds, type="TPM")
lapply(names(markers_exp_list), function(name) {
  file_name <- file.path(fig_dir, paste0("Year2_Boxplot_Markers_", name, ".pdf"))
  data <- markers_exp_list[[name]] %>%
    mutate(log.value=log10(value+1))
  gg <- ggplot(data, aes(x=group, y=log.value, group=group)) +
    geom_boxplot(outlier.shape = NA) +
    #geom_violin() +
    geom_jitter(width=0.2, size=1.2) +
    # geom_point() +
    facet_wrap( ~ gene_name, nrow=1, scale="free") +
    theme_bw() +
    #scale_y_log10() +
    labs(y=bquote(log[10]~"(TPM+1)")) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
          axis.title.x = element_blank(),
          legend.justification=c(0,1), legend.position=c(0, 1)) 

  pdf(file_name, height=3, width=7)
  plot(gg)
  dev.off()         
})


#'
#' relationship with MRI/Pathology and DUX4 scores

#'
#' (4) correlation (Spearman) of RNA-seq to MRI and histopath changes 
#'  use corrgram and update to figures that uses in previous paper
dux4_score <- get(load(file.path(pkg_dir, "public_data", "year2_dux4_score.rda")))
time_2 <- mri_pathology[["time_2"]] %>% 
  left_join(as.data.frame(dux4_score), by="sample_name") %>% 
  dplyr::filter(!is.na(Pathology.Score)) %>%
  rename(DUX4.Score = dux4.rlogsum)
    

#' (4b) scatter of DUX4 and pathology score
cor <- cor(time_2$Pathology.Score, time_2$DUX4.Score,
           method="spearman")
gg <- ggplot(time_2, aes(x=Pathology.Score, y=DUX4.Score))+
  geom_point(size=1) +
  geom_smooth(method = lm, se=TRUE, color="grey70", linetype="dashed") + 
  annotate("text", x=7.5, y=15, vjust=0.5, 
           label=paste0("correlation=", format(cor, digits=3))) +
  theme_bw() +
  labs(y="DUX4 score \n (relative biomarker expression)",
       x="Pathology score") 
pdf(file.path(fig_dir, "Year2_DUX4score_Patho_ScatterPlot.pdf"), 
    width=4, height=3)
plot(gg)
dev.off()


#'
#' MRI vs DUX4.Score
#'
dux4_score <- get(load(file.path(pkg_dir, "public_data", "year2_dux4_score.rda")))
time_2 <- mri_pathology[["time_2"]] %>% 
  left_join(as.data.frame(dux4_score), by="sample_name") %>% 
  rename(DUX4.Score = dux4.rlogsum) 

#' (4d) Scatter: DUX4.score vs STIR rating
cor <- cor(as.numeric(time_2$STIR_rating), time_2$DUX4.Score,
           method="spearman")
gg <- ggplot(time_2, aes(x=as.numeric(STIR_rating), y=DUX4.Score))+
  geom_point(size=1) +
  geom_smooth(method = lm, se=TRUE, linetype="dashed", color="grey70") + 
  annotate("text", x=2.5, y=14.5, vjust=0.5, 
           label=paste0("correlation=", format(cor, digits=3))) +
  theme_bw() +
  labs(y="DUX4 score \n (relative biomarker expression)",
       x="STIR rating") 
pdf(file.path(fig_dir, "Year2_DUX4score_STIR_ScatterPlot.pdf"), 
    width=4, height=3)
plot(gg)
dev.off()

#' (4d) Boxplot: DUX4.score vs STIR+/MRI normal 
#' Note: new_time_2 is made for STIR+/MRI normal plot
#' 1614b and 32-00007b with MRI normal but high DUX4.score
DUX4_threshold <- 2.5
new_time_2 <- time_2 %>%
  mutate(STIR_status = ifelse(STIR_rating > 0, "STIR+", "MRI normal")) %>%
  mutate(DUX4_status = ifelse(DUX4.Score > DUX4_threshold, "DUX4+", "DUX4-")) %>%
  mutate(STIR_status = ifelse((STIR_rating == 0 & T1_rating > 0), NA, STIR_status)) %>%
  dplyr::filter(!is.na(STIR_status))

gg_stir <- ggplot(new_time_2, aes(x=STIR_status, y=DUX4.Score,
                       group=STIR_status)) +
  geom_boxplot(outlier.shape = NA, width=0.5 ) +
  geom_jitter(width=0.2, size=1) +
  geom_hline(yintercept = DUX4_threshold, color="grey", 
             linetype="dashed") +
  theme_bw() + 
  labs(title="MRI normal vs. STIR+",
       y="DUX4 score", x="STIR status")  +
  theme(#axis.title.x = element_blank(),
        #axis.title.y = element_text(size=9),
        plot.title = element_text(size=9, hjust = 0.5))
pdf(file.path(fig_dir, "Year2_DUX4score_STIR_Boxplot.pdf"),
    width=3, height=3)  
print(gg_stir)
dev.off()


#' 1614b and 32-00007b with MRI normal but high DUX4.score
#' some statistics calculation
new_time_2 %>% filter(STIR_status == "MRI Normal")
new_time_2 %>% group_by(STIR_status) %>% summarize(avg=mean(DUX4.Score), sd=sd(DUX4.Score))
new_time_2 %>% filter(STIR_status == "STIR+") %>%
  summarize(pert = sum(DUX4_status == "DUX4+") / 18) #total 18 rows of STIR+ (72.2% DUX4+ )

#' (4f) DUX4 score and T1 fraction (fat_fraction) 

gg <- ggplot(time_2, aes(x=fat_fraction, y=DUX4.Score))+
  geom_point(size=1) +
  geom_smooth(method = loess, se=FALSE, linetype="dashed", color="gray50") + 
  geom_vline(xintercept = 0.2, linetype="dashed", color="gray") +
  geom_vline(xintercept = 0.4, linetype="dashed", color="gray") +
  theme_bw() +
  labs(y="DUX4 score \n (relative biomarker expression)",
       x="T1 fat fraction")     
pdf(file.path(fig_dir, "Year2_DUX4score_FF_ScatterPlot.pdf"), 
    width=4, height=3)
plot(gg)
dev.off()

cor <- cor(time_2$fat_fraction, time_2$DUX4.Score, method="spearman")
gg <- ggplot(time_2, aes(x=fat_fraction, y=DUX4.Score))+
  geom_point(size=1) +
  geom_smooth(method = lm, se=TRUE, linetype="dashed", color="gray70") + 
  annotate("text", x=0.5, y=10, vjust=0.5, 
           label=paste0("correlation=", format(cor, digits=3))) +
  theme_bw() +
  labs(y="DUX4 score \n (relative biomarker expression)",
       x="T1 fat fraction")     
pdf(file.path(fig_dir, "Year2_DUX4score_FF_ScatterPlot_linear.pdf"), 
    width=4, height=3)
plot(gg)
dev.off()


new_time_2 <-  time_2 %>%
  mutate(T1_fraction_bins = cut(fat_fraction, breaks=c(0, 0.2, 0.25, 0.4, 0.5, 0.72)))
gg <- ggplot(new_time_2, aes(x=T1_fraction_bins, y=DUX4.Score)) +
  geom_boxplot(outlier.shape = NA, width=0.5) +
  #geom_point(size=0.8) +
  labs(x="T1 fraction bins", y="DUX4 score") +
  geom_jitter(width=0.2, size=0.8) +
  theme_bw() +
  theme(#axis.title.x = element_blank(),
        #axis.title.y = element_text(size=9),
        plot.title = element_text(size=9, hjust = 0.5))
pdf(file.path(fig_dir, "Year2_DUX4score_FF_bins_boxPlot.pdf"), 
    width=4, height=3)
plot(gg)
dev.off()

#' inflammation infiltrates vs DUX4 score
gg <- ggplot(time_2, aes(x=Inflammation, y=DUX4.Score)) +
  geom_boxplot(outlier.shape=NA, width=0.5) +
  geom_jitter(width=0.2, size=0.8) +
  theme_bw()
pdf(file.path(fig_dir, "Year2_DUX4score_inflamm.pdf"), width=4, height=3)
plot(gg)
dev.off()  



#'############################
#' supplemental figures
#'###########################

#'
#' (A) loading variables of PC1 and PC2 and associated GO terms
#'
pkg_dir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsy.all"
load(file.path(pkg_dir, "public_data", "year2.dds.rda"))
load(file.path(pkg_dir, "public_data", "year2.rlg.rda"))
library(DESeq2)

#' get PCA
pca <- getPCA(year2.rlg, ntop=500)
pcs <- pca$x
loading_var <- as.data.frame(pca$rotation)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
names(percentVar) <- colnames(pca$x)

#' get markers
library(xlsx)
high.class <- 
  read.xlsx(file=file.path(pkg_dir, "manuscript", "tables", 
                           "Suppl_table_6_Candidate_Biomarkers_and_Enriched_Go.xlsx"),
            sheetIndex=1)
get_de_marker <- function(df, col_name) {
  col_name <- sym(col_name)
  df %>% dplyr::filter(!!col_name==TRUE) %>% pull(gencode_id) %>% as.character()
}            
marker_flags <- list(
  dux4 = get_de_marker(high.class, "DUX4_induced"),
  inflamm = get_de_marker(high.class, "inflamm"),
  immune = get_de_marker(high.class, "immune"),
  ecm  = get_de_marker(high.class, "ecm"))

#' get loading variables and tidy the data.frame
loading_var <- as.data.frame(pca$rotation) %>%
  rownames_to_column(var="gencode_id") %>%
  dplyr::select(gencode_id, PC1, PC2) %>%
  mutate(gene_name = rowData(year2.rlg[as.character(gencode_id)])$gene_name) %>%
  mutate(functions=ifelse(gencode_id %in% marker_flags$inflamm, "inflam", NA)) %>%
  mutate(functions=ifelse( is.na(functions) & (gencode_id %in% marker_flags$immune), "immune", functions)) %>%
  mutate(functions=ifelse( is.na(functions) & (gencode_id %in% marker_flags$ecm), "ecm", functions)) %>%
  mutate(functions=ifelse( is.na(functions) & (gencode_id %in% marker_flags$dux4), "DUX4", functions)) 
  
pc1_loading <- loading_var %>%
  arrange(desc(abs(PC1))) %>%
  mutate(gene_name = factor(gene_name, levels=as.character(gene_name))) 
tmp = ifelse(grepl("IGH|IGK|IGL", pc1_loading$gene_name) & is.na(pc1_loading$functions), "Ig", pc1_loading$functions)
tmp = ifelse(grepl("MYH|MYL", pc1_loading$gene_name) & is.na(tmp), "muscle", tmp)
pc1_loading$functions <- tmp

pc2_loading <- loading_var %>%
  arrange(desc(abs(PC2))) %>%
  mutate(gene_name = factor(gene_name, levels=as.character(gene_name)))
tmp =ifelse(grepl("IGH|IGK|IGL", pc2_loading$gene_name) & is.na(pc2_loading$functions), "Ig", pc2_loading$functions)
tmp = ifelse(grepl("MYH|MYL", pc2_loading$gene_name) & is.na(tmp), "muscle", tmp)
pc2_loading$functions <- tmp

#' barplots
gg_pc1 <- ggplot(pc1_loading[1:60, ], aes(x=gene_name, y=PC1, fill=functions)) +
  geom_bar(stat="identity") +
  labs(title="Top 60 PC1 loading variable", y="rotation") +
  theme_bw() + 
  theme(legend.key.size = unit(0.8, "line"), 
    axis.title.x=element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1,
                                   size=6))
gg_pc2 <- ggplot(pc2_loading[1:60, ], aes(x=gene_name, y=PC2, fill=functions)) +
  geom_bar(stat="identity") +
  labs(title="Top 60 PC2 loading variable", y="rotation") +
  theme_bw() + 
  theme(legend.key.size = unit(0.8, "line"),
    axis.title.x=element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1,
                                   size=6))                                   
library(gridExtra)
pdf(file.path(sub_dir, "year2_pca_loading_vars.pdf"), width=8, height=6)
grid.arrange(gg_pc1, gg_pc2, nrow=2)
dev.off()

#'
#' (B) PCA vis for all year II and show that 32-0008b is an outlier
#'
sub.dds <- sanitized.dds[, 
   (sanitized.dds$visit == "II" | sanitized.dds$pheno_type == "Control")]
#' re-normalized again
sub.dds <- estimateSizeFactors(sub.dds)
sub.dds <- estimateDispersions(sub.dds)
design(sub.dds) <- ~ pheno_type
sub.rlg <- rlog(sub.dds, blind = FALSE) #normalized with respect to year2 library size

data <- plotPCA(sanitized.rlg, 
                intgroup=c("pheno_type", "visit"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
gg <- ggplot(data, aes(PC1, PC2, color=pheno_type, shape=visit)) +
  geom_point(size=1.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  labs(title="PCA of both years including outlier 32-0008b") +
  #geom_text_repel(aes(label=rownames(data)), size=3,  show.legend=FALSE) +
  geom_text(x=50, y=85, label="32-0008b", hjust=0, 
           show.legend=FALSE, size=3) +
  theme_minimal() +
  theme(legend.position="bottom", 
        plot.title = element_text(hjust = 0.5, size=10))
pdf(file.path(sup_dir, "SuppFigA_Year2_PCA.pdf"), width=6, height=4)
plot(gg)
dev.off()  

#' (C) Heatmap of 53 biomarkers on group 1 to 4 (zscore of rlog)
library(pheatmap)
dux4_targets <- get(load(file.path(pkg_dir, "data", "FSHD_markers.rda")))
dux4_targets <- as.data.frame(dux4_targets) %>% 
  filter(gene_name != "TRIM43CP")
year2.dds$dux4.group <- factor(year2.dds$dux4.group)
data_lst <- lapply(levels(year2.dds$dux4.group), function(x, rlg) {
  s <- colnames(rlg)[rlg$dux4.group == x]
  data <- assay(rlg[dux4_targets$gene_id, s]) 
  # re-arrange by the colmeans (look heatmap)
  Colv <- colMeans(data, na.rm = TRUE)
  data <- data[, order(Colv)]
}, rlg=year2.rlg)
data <- do.call(cbind, data_lst)
rownames(data) <- rowData(year2.dds[rownames(data)])$gene_name
zscore_data <- (data - rowMeans(data)) / rowSds(data)

annotation_col <- 
  data.frame(dux4.group = year2.dds[, colnames(data)]$dux4.group,
             phenotype = year2.dds[, colnames(data)]$pheno_type)
rownames(annotation_col) <- colnames(data)
gaps_col <- cumsum(count(annotation_col, dux4.group)[, "n"])$n

pheatmap(data, annotation_col=annotation_col,
         fontsize_row=7, fontsize_col=6,
         scale="none", 
         gaps_col=gaps_col,
         #annotation_colors=ann_cor,
         cellheight=8, #annotation_legend=FALSE,
         cluster_rows=TRUE, cluster_cols=FALSE,
         file=file.path(sup_dir, "year2_53markers_heatmap.pdf"))
