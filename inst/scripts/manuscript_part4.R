#' this script 
#' (1) exam mild class vs controls and gives biomarkers for detecting early muscle changes
#' (2) PAX7 and PAX3 expression
#'

#' load library and datasets
library(DESeq2)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(BiocParallel)
multi_param <- MulticoreParam(worker=4)
register(multi_param, default=TRUE)

pkg_dir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsy.all"
fig_dir <- file.path(pkg_dir, "manuscript", "figures")
sub_dir <- file.path(pkg_dir, "manuscript", "sup_figs")
table_dir <- file.path(pkg_dir, "manuscript", "tables")
setwd(pkg_dir)

source(file.path(pkg_dir, "scripts", "tools.R"))
source(file.path(pkg_dir, "scripts", "manuscript_tools.R"))
load(file.path(pkg_dir, "public_data", "cluster_df.rda"))
load(file.path(pkg_dir, "public_data", "sanitized.dds.rda"))
load(file.path(pkg_dir, "public_data", "sanitized.rlg.rda"))

sanitized.dds$RNA_cluster_5 <- sanitized.rlg$RNA_cluster_5 <- cluster_df$RNA_cluster
sanitized.dds$new_cluster_name <- sanitized.rlg$new_cluster_name <- cluster_df$new_cluster_name

mk <- get(load(file.path(pkg_dir, "data", "FSHD_markers.rda")))$gene_name
mk <- mk[mk %in% rowData(sanitized.dds)$gene_name]

.tidy_results <- function(res, dds, padj_thres=0.05) {
  de <- as.data.frame(res) %>%
    rownames_to_column(var="gencode_id") %>%
    dplyr::filter(padj < padj_thres) %>%
    mutate(gene_name = rowData(dds[gencode_id])$gene_name) %>%
    mutate(gene_id =
             sapply(strsplit(gencode_id, ".", fixed=TRUE), "[[", 1))
}

#'
#' (1) differential analysis for Mild (A) vs A_Cntr
#'
new.dds <- sanitized.dds
design(new.dds) <- ~ RNA_cluster_5
new.dds <- DESeq(new.dds)
deseq2_design_cluster <- new.dds
save(deseq2_design_cluster, file=file.path(pkg_dir, "data", "deseq2_design_cluster.rda"))
save(deseq2_design_cluster, file=file.path(pkg_dir, "public_data", "deseq2_design_cluster.rda"))

#'
#' (2) Get Differentially expressed genes in class A (Mild) vs. A_Cntr
#'
new.dds <- get(load(file.path(pkg_dir, "public_data", "deseq2_design_cluster.rda")))
classA.res <- results(new.dds, alpha = 0.05,
                      name = "RNA_cluster_5_A_vs_A_Cntr")
de_A <- .tidy_results(classA.res, dds=new.dds, padj=0.05)                      
de_A_up <- de_A %>% dplyr::filter(log2FoldChange > 0) %>%
  mutate(DUX4_marker=gene_name %in% mk)
de_A_down <- as.data.frame(classA.res) %>% 
  rownames_to_column(var="gene_id") %>%
  filter(padj < 0.05, log2FoldChange < 0) %>%
  mutate(gene_name=rowData(sanitized.dds[gene_id])$gene_name)
#: class D
classD.res <- results(new.dds, alpha = 0.05, lfcThreshold=2,
                      name="RNA_cluster_5_D_vs_A_Cntr")          
de_D_up <- .tidy_results(classD.res, dds=new.dds, padj=0.05) %>%
   dplyr::filter(log2FoldChange > 0) %>%
  mutate(DUX4_marker=gene_name %in% mk)

#'
#' Heatmap of de_A_up 
#'
data <- .get_local_rlg(de_A_up$gencode_id, sanitized.dds) 
rownames(data) <- rowData(sanitized.dds[rownames(data)])$gene_name
#rownames(data) <- de_A_up$gene_name
ann_col <- data.frame(pheno_type = new.dds$pheno_type,
                      cluster = new.dds$RNA_cluster_5)
rownames(ann_col) <- colnames(data)
zscore_data <- (data - rowMeans(data)) / rowSds(data)
pheatmap(zscore_data, annotation_col=ann_col,
         fontsize_row=5, 
         #fontsize_col=6,
         scale="none", 
         show_rownames=TRUE,
         cellheight=5, 
         file=file.path(sub_dir, "mild_FSHD_up_regulated.pdf"))

#' only A and A_Cntr
idx <- sanitized.dds$RNA_cluster_5 %in% c("A", "A_Cntr")
data <- .get_local_rlg(de_A_up$gencode_id, sanitized.dds[, idx])
rownames(data) <- rowData(sanitized.dds[rownames(data)])$gene_name
#rownames(data) <- de_A_up$gene_name
ann_col <- data.frame(pheno_type = new.dds$pheno_type[idx])
rownames(ann_col) <- colnames(data)
zscore_data <- (data - rowMeans(data)) / rowSds(data)
pheatmap(zscore_data, annotation_col=ann_col,
         fontsize_row=5, 
         #fontsize_col=6,
         scale="none", 
         show_rownames=TRUE,
         cellheight=5, 
         file=file.path(sub_dir, "mild_FSHD_up_regulated_clusterA.pdf"))

#' 
#' (1) select de_up that are robust in cluster D (High, lfc > 2 with padj < 0.05)
#' (2) define key predictors: limited to DE genes in High (class D)
#' (3) plot heatmap to see the progression of expression in other FSHD classes
#'
keys <- de_A_up %>% dplyr::filter(gencode_id %in% de_D_up$gencode_id) %>%
  arrange(padj)
  
data <- assay(sanitized.rlg[keys$gencode_id, ])
rownames(data) <- keys$gene_name
ann_col <- data.frame(pheno_type = new.dds$pheno_type,
                      cluster = new.dds$RNA_cluster_5)
rownames(ann_col) <- colnames(data)
zscore_data <- (data - rowMeans(data)) / rowSds(data)

pheatmap(zscore_data, annotation_col=ann_col,
         fontsize_row=5, 
         fontsize_col=6,
         scale="none", 
         show_rownames=TRUE,
         cellheight=5, 
         file=file.path(sub_dir, "mild_164_key_predictors_allsamples.pdf"))

#' run goseq on these 164 genes


#'
#' goseq for key predictors: lymphocyte activation, T/B cell...
#
keys <- de_A_up %>% dplyr::filter(gencode_id %in% de_D_up$gencode_id) %>%
  arrange(padj)
keep <- sanitized.rlg$RNA_cluster_5 %in% c("A_Cntr", "A")
tmp.rlg <- sanitized.rlg[keys$gencode_id, keep]
tmp.dds <- sanitized.dds[keys$gencode_id, keep]
tmp.rlg$RNA_cluster_5 <- factor(tmp.rlg$RNA_cluster_5, levels=c("A_Cntr", "A"))
rownames(tmp.rlg) <- rowData(tmp.rlg)$gene_name

#' goseq for the keys
universal <- sapply(strsplit(rownames(classA.res), ".", fixed=TRUE),
                    "[[", 1)
enriched <- .do_goseq(universal = universal, 
                      selected_genes = keys$gene_id,
                      return.DEInCat=FALSE, dds=sanitized.dds)

#' (1) there are two groups of Class A: 2422, 32-0019, 1614, 32-0018, 32-0017b (Mild_DUX4+) have elevated
#' expression in PRAMEFX (DUX4-induced genes) whereas the other group (Mild_DUX4-) do not have those
#' expression. Note 32-0014 is highly elevated in some genes. 

library(randomForest)
library(MLInterfaces)
library(ROC)
library(genefilter)
set.seed(123)
data <- .get_local_rlg(keys$gencode_id, tmp.dds)
rownames(data) <- rowData(tmp.dds[rownames(data)])$gene_name
#rownames(data) <- keys$gene_name
exp_set <- ExpressionSet(data,
                         phenoData   = AnnotatedDataFrame(as.data.frame(colData(tmp.dds))))

#' pick training set: use follow-up samples as test sets (let's do leave-one-out )
train_ind <- c(1:11, 14:19)
rf1 <- MLearn(pheno_type ~., data=exp_set,
              randomForestI, train_ind, ntree=1000, mtry=10, importance=TRUE) 
confuMat(rf1, "train")              
confuMat(rf1, "test")

#' leave-one-out cross-validation
error_test <- lapply(1:ncol(exp_set), function(i) {
  train_ind <- c(1:25)[-i]
  rf <- MLearn(pheno_type ~., data=exp_set,
              randomForestI, train_ind, ntree=1000, mtry=10, importance=TRUE) 
  cfm <- confuMat(rf, "test")
  error_sample <- ifelse(is.table(cfm), 1, 0)
})
error_rate <- sum(unlist(error_test)) / ncol(exp_set)

#' prediction error rate is almost 0
imp <- getVarImp(rf1) # rownames are messed up by MLearn
rownames(imp) <- featureNames(exp_set)

tmp_data <- report(imp, n=25) %>% 
  arrange(MnDecrAcc) %>%
  mutate(names=as.character(names)) %>%
  mutate(names=factor(names, levels=names))

gg_imp <- ggplot(tmp_data, aes(x=names, y=MnDecrAcc)) +
  geom_bar(stat="identity") + 
  coord_flip() +
  theme_bw()

imp_data <- as.data.frame(imp@.Data) %>%
  rownames_to_column(var="gene_name") %>%
  #left_join(keys, by="gene_name") %>%
  arrange(desc(MeanDecreaseAccuracy)) 

#'
#' ROC
#'
rocs <- rowpAUCs(exp_set, "pheno_type", p=0.2)
rocs_data <- data.frame(pAUC = area(rocs), AUC = rocs@AUC) %>%
  rownames_to_column(var="gene_name")                        
                    
comb_data <- left_join(imp_data, rocs_data, by="gene_name") %>% 
  #dplyr::filter(pAUC > 0.1) %>%
  dplyr::select(-Control, -FSHD, -MeanDecreaseGini) %>%
  left_join(keys %>% dplyr::select(-c(gencode_id, gene_id)), by="gene_name")


write.csv(comb_data, file=file.path(table_dir, "mild_FSHD_potential_markers.csv"),
          row.names=FALSE)

#'
#' cross-validation using AUC > 0.9
#' 
sub_exp <- exp_set[rocs_data %>% dplyr::filter(AUC>0.9) %>% pull(gene_name)]
error_test <- lapply(1:ncol(sub_exp), function(i) {
  train_ind <- c(1:25)[-i]
  rf <- MLearn(pheno_type ~., data=sub_exp,
              randomForestI, train_ind, ntree=1000, mtry=10, importance=TRUE) 
  cfm <- confuMat(rf, "test")
  error_sample <- ifelse(is.table(cfm), 1, 0)
})

error_rate <- sum(unlist(error_test)) / ncol(sub_exp)

#' heatmap
predictors <- comb_data %>% dplyr::filter(AUC >= 0.9) %>% 
  mutate(gencode_id = get_ensembl(gene_name, sanitized.dds)) %>%
 pull(gencode_id) #52

keep <- sanitized.rlg$RNA_cluster_5 %in% c("A_Cntr", "A")
tmp.rlg <- sanitized.rlg[predictors, keep]
tmp.dds <- sanitized.dds[predictors, keep]
#data <- assay(tmp.rlg)
data <- .get_local_rlg(predictors, tmp.dds)
rownames(data) <- rowData(tmp.dds)$gene_name
ann_col <- data.frame(pheno_type=tmp.dds$pheno_type,
                      class = tmp.dds$RNA_cluster_5)
rownames(ann_col) <- colnames(data)
zscore_data <- (data - rowMeans(data)) / rowSds(data)
pheatmap(zscore_data, annotation_col=ann_col,
         fontsize_row=6.5, 
         #fontsize_col=6,
         scale="row", 
         show_rownames=TRUE,
         cellheight=8, 
         file=file.path(fig_dir, "mild_FSHD_key_predictors_heatmap.pdf"))


#'
#' visualization for some cadindates
#' pick top 25 importance and pAUC > 0.1
#'
imp_data %>% top_n(25) %>% filter(pAUC > 0.15) %>% arrange(pAUC)

sel <- c("CDKN1A", "C1QA", "CD4", "TNC", "RUNX1", "LILRB5", "FGF18", "CCL18")
pdf(file.path(sub_dir, "mild_FSHD_supervised_ML.pdf"))  
plot(gg_imp)
lapply(sel, function(name) plot(rocs[name], main=name))
dev.off()

#' boxplot of TPM 
factor <- cluster_df$new_name
names(factor) <- as.character(cluster_df$sample_name)
set_exp <- .getExpFromDDS_by_cluster(sel, dds=sanitized.dds, factor=factor, type="TPM")
data <- set_exp %>%
  mutate(log.value=log10(value+1))
gg <- ggplot(data, aes(x=group, y=log.value, group=group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width=0.2, size=0.7) +
    facet_wrap( ~ gene_name, nrow=2, scale="free") +
    theme_bw() +
    labs(y=bquote(log[10]~"(TPM+1)")) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
          axis.title.x = element_blank(),
          legend.justification=c(0,1), legend.position=c(0, 1)) 

pdf(file.path(sub_dir, "mild_FSHD_some_predictor.pdf"))
plot(gg)
dev.off() 

#plot(rocs[name], main=name))   

#' CDKN1A ROC plot and boxplot
gg_CDKN1A <- ggplot(data %>% filter(gene_name=="CDKN1A"), aes(x=group, y=log.value, group=group)) +
    geom_boxplot(outlier.shape = NA, width=0.7) +
    #geom_violin() +
    geom_jitter(width=0.2, size=0.8) + 
    labs(y=bquote(log[10]~"(TPM+1)"), title="CDKN1A") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust=0.5, size=10),
          axis.title.y = element_text(size=9),
          legend.justification=c(0,1), legend.position=c(0, 1)) 
pdf(file.path(fig_dir, "mild_FSHD_CDKN1A_vs_classes.pdf"),
    width=3.5, height=3.5)
plot(gg_CDKN1A)
dev.off()  
pdf(file.path(fig_dir, "mild_FSHD_CDKN1A_ROC.pdf"), width=6, height=6)
plot(rocs["CDKN1A"], main="CDKN1A")
dev.off()

#' CCNA1 ROC plot and boxplot
factor <- as.character(cluster_df$new_cluster_name)
factor[cluster_df$RNA_cluster=="A_Cntr"] <- "Control"
factor <- factor(factor, levels=c("Control", "Mild", "Moderate", "IG-High", "High", "Muscle-Low"))
names(factor) <- as.character(cluster_df$sample_name)
ccna1_exp <- .getExpFromDDS_by_cluster("CCNA1", dds=sanitized.dds, factor=factor, type="TPM")
data <- ccna1_exp %>%
  mutate(log.value=log10(value+1))
gg_CCNA1 <- ggplot(data %>% filter(gene_name=="CCNA1"), aes(x=group, y=log.value, group=group)) +
    geom_boxplot(outlier.shape = NA, width=0.7) +
    #geom_violin() +
    geom_jitter(width=0.2, size=0.8) + 
    labs(y=bquote(log[10]~"(TPM+1)"), title="CCNA1") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust=0.5, size=10),
          axis.title.y = element_text(size=9),
          legend.justification=c(0,1), legend.position=c(0, 1))   
pdf(file.path(fig_dir, "mild_FSHD_CCNA1_vs_classes.pdf"),
    width=3.5, height=3.5)
plot(gg_CCNA1)
dev.off()  
pdf(file.path(fig_dir, "mild_FSHD_CCNA1_ROC.pdf"), width=6, height=6)
plot(rocs["CCNA1"], main="CCNA1")
dev.off()
#'
#' predictor's value over two points


#'
#' how about PAX7 and PAX3
#'

#' What's PAX3 and PAX7  in each of the classes? plot boxplot
#' What's PAX3 and PAX7 ROC and importance by Random forest???  not what I expected
#' what's logFC of PAX4in Mild FSHD and Hot vs Control? -0.443319961787262, -1.00196639288724
#' what's logFC of PAX7 in Mild FSHD and Hot vs Control? -0.223548152412414, 0.458500515433963

pax_id <- get_ensembl(c("PAX3", "PAX7"), sanitized.dds)      
classD.res[pax_id, ]
classA.res[pax_id, ]

tmp.dds <- sanitized.dds[, keep]
data <- .get_local_rlg(c(pax_id, keys$gencode_id), tmp.dds)
rownames(data) <- rowData(tmp.dds[c(pax_id, keys$gencode_id)])$gene_name
exp_set <- ExpressionSet(data,
                         phenoData = AnnotatedDataFrame(as.data.frame(colData(tmp.dds))))
#' pick training set: use follow-up samples as test sets
train_ind <- c(1:12, 17:22)
rf1 <- MLearn(pheno_type ~., data=exp_set,
              randomForestI, train_ind, ntree=1000, mtry=10, importance=TRUE) 
confuMat(rf1, "train")              
confuMat(rf1, "test")
imp <- getVarImp(rf1)     
as.data.frame(imp@.Data) %>%
  rownames_to_column(var="gene_name") %>%
  arrange(MeanDecreaseAccuracy) # PAX7 is at #7

#' PAX7's ROC?  
rocs <- rowpAUCs(exp_set, "pheno_type", p=0.2)
rocs_data <- data.frame(pAUC = area(rocs), AUC = rocs@AUC) %>%
  rownames_to_column(var="gene_name") %>%
  arrange(desc(pAUC)) # PAX7 is #49, as bad as KHDC1L, PAX3 is #37
rocs_data

#' boxplot with selected predictors, including PAX7
sel <- c("PAX7", "PAX3")
factor <- cluster_df$new_name
names(factor) <- as.character(cluster_df$sample_name)
set_exp <- .getExpFromDDS_by_cluster(sel, dds=sanitized.dds, factor=factor, type="TPM")
data <- set_exp %>%
  mutate(log.value=log10(value+1))
gg <- ggplot(data, aes(x=group, y=log.value, group=group)) +
    geom_boxplot(outlier.shape = NA) +
    #geom_violin() +
    geom_jitter(width=0.2, size=0.7) +
    # geom_point() +
    facet_wrap( ~ gene_name, nrow=1, scale="free") +
    theme_bw() +
    #scale_y_log10() +
    labs(y=bquote(log[10]~"(TPM+1)")) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=10),
          legend.justification=c(0,1), legend.position=c(0, 1)) 

pdf(file.path(sub_dir, "PAX7_boxplot_by_cluster.pdf"), width=4, height=3.5)
plot(gg)
dev.off()         


#'
#' more abourt PAX7: PAX7-supressed and DUX4-induced
#'
library(xlsx)
sheet_file <- file.path(pkg_dir, "extdata", "PAX7", "41467_2017_1200_MOESM6_ESM.xlsx")
up <- read.xlsx(sheet_file, sheetIndex=1, startRow=3, colIndex=c(1:3), stringsAsFactors=FALSE) %>%
  add_column(attribute="up_regulated") 

down <- read.xlsx(sheet_file, sheetIndex=1, startRow=3, colIndex=c(4:6), stringsAsFactors=FALSE) %>%
  add_column(attribute="down_regulated")

i <- de_A_up$gene_name %in% $Human.Gene.Symbol 
#' heatmap of down$Human.Gene.Symbol
de_A_up %>% filter(gene_name %in% down$Human.Gene.Symbol)
tmp <- intersect(down$Human.Gene.Symbol, rowData(sanitized.dds)$gene_name)
pax7_supressed_id <- get_ensembl(tmp, sanitized.dds)

tmp <- intersect(up$Human.Gene.Symbol, rowData(sanitized.dds)$gene_name)
pax7_induced_id <- get_ensembl(tmp, sanitized.dds)
# supressed vs mild fshd
as.data.frame(classA.res[pax7_supressed_id, ]) %>%
  rownames_to_column(var="gencode_id") %>%
  #summarize(mean(log2FoldChange))
  filter(log2FoldChange > 1) %>%
  mutate(gene_name = rowData(sanitized.dds[gencode_id])$gene_name) %>%
  pull(gene_name)

de_D_up %>% filter(gene_name %in% down$Human.Gene.Symbol) 
as.data.frame(classD.res[pax7_supressed_id, ]) %>%
  rownames_to_column(var="gencode_id") %>%
  #summarize(mean(log2FoldChange))
  filter(log2FoldChange > 1) %>%
  mutate(gene_name = rowData(sanitized.dds[gencode_id])$gene_name)  
  


# induced vs mild
de_A_up %>% filter(gene_name %in% up$Human.Gene.Symbol) %>% # SLC37A2
# eight have log2FC > 1 but will poor adjusted p-value
as.data.frame(classA.res[pax7_induced_id, ]) %>%
  rownames_to_column(var="gencode_id") %>%
  summarize(mean(log2FoldChange))
  filter(log2FoldChange > 1) %>%
  mutate(gene_name = rowData(sanitized.dds[gencode_id])$gene_name) 

de_D_up %>% filter(gene_name %in% up$Human.Gene.Symbol)
as.data.frame(classD.res[pax7_induced_id, ]) %>%
  rownames_to_column(var="gencode_id") %>%
  #summarize(mean(log2FoldChange))
  filter(log2FoldChange > 1) %>%
  mutate(gene_name = rowData(sanitized.dds[gencode_id])$gene_name) 
