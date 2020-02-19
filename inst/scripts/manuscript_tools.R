get_ensembl <- function(inquiry_genename, rse) {
  #' convert inquiry_genename to gene_id using embedding annotation
  #' in rse 
  #' Note: rse is the summarizedExperiment-type of instance that
  #' obtain annotation in rowData
  require(dplyr)
  #' sanity check
  inquiry_genename <- unique(inquiry_genename)
  flag <- all(inquiry_genename %in% rowData(rse)$gene_name)
  #cat(inquiry_genename, "\n")
  if (!flag)
    warning("Not all inquiry matched the embedded annotation in rse")
  
  if (flag) {
    id <- as.data.frame(rowData(rse)) %>% 
      dplyr::select(gene_name, gene_id) %>%
      dplyr::filter(gene_name %in% inquiry_genename) %>%
      column_to_rownames(var="gene_name")
    id <- id[inquiry_genename, ]
    names(id) <- inquiry_genename
    return(id)
  }
  
}

.get_DUX4_clust_tree <- function(dux4_rlog) {
    # for year2 data
    #dux4_rlog <- assay(rlg)
    dux4_hc <- hclust(dist(t(dux4_rlog)))
    dux4_tree <- cutree(dux4_hc, k=3)

    #' 5> dux4.score > 1.5 samples to 2 (01-0022b, 01-0030b and 32-0009b)
    # boost 2 -> 3, 3 -> 4
    group1_idx <- which(dux4_tree == 1)
    dux4_tree <- dux4_tree + 1
    dux4_tree[group1_idx] <- 1
    dux4_tree[c("01-0022b", "01-0030b", "32-0009b")] <- 2
    dux4_tree
}

.insert_gene_name <- function(exp, sub_dds) {
  exp <- as.data.frame(exp) %>%
    rownames_to_column(var="gencode_id") 
  rownames(exp) <- rowData(sub_dds[exp$gencode_id])$gene_name
  exp
}

.write_expression_to_xlsx <- function(sub_dds, sub_rlg=NULL,
                                      file_name,
                                      append=TRUE) {
  #' sub_dds is DESeq2's dds instance, append genename as rownames
  cnt <- counts(sub_dds, normalized=FALSE)
  cnt <- .insert_gene_name(cnt, sub_dds=sub_dds)

  norm.cnt <- counts(sub_dds, normalized=TRUE)
  norm.cnt <- .insert_gene_name(norm.cnt, sub_dds=sub_dds)

  tpm <- assays(sub_dds)$TPM
  tpm <- .insert_gene_name(tpm, sub_dds=sub_dds)

  rpkm <- assays(sub_dds)$RPKM
  rpkm <- .insert_gene_name(rpkm, sub_dds=sub_dds)

  write.xlsx(cnt, sheetName="counts", 
             file=file_name, append=TRUE)
  write.xlsx(norm.cnt, sheetName="normalized counts",
             file=file_name, append=TRUE)
  write.xlsx(tpm, sheetName="TPM",
             file=file_name, append=TRUE)
  write.xlsx(rpkm, sheetName="RPKM",
             file=file_name, append=TRUE)            
}

.makeFourBiomarkerGroupTable <- function(dds, relative_logsum, 
                                         dux4_group,
                                         markers_id, file_name) {
  sub_dds <- dds[markers_id]
  sample_name <- names(relative_logsum)

  output <- as.data.frame(relative_logsum) %>%
    dplyr::rename(dux4.rlogsum=relative_logsum) %>%
    rownames_to_column(var="sample_name") %>% 
    mutate(dux4.score=dux4.rlogsum - min(dux4.rlogsum)) %>% 
    add_column(dux4.group=dux4_group[sample_name]) %>% 
    add_column(pheno_type=sub_dds[ ,sample_name]$pheno_type,
               .before="dux4.rlogsum")

  #' get count data
  cnt <- as.data.frame(t(counts(sub_dds, normalized=FALSE)))
  colnames(cnt) <- names(markers_id)
  norm.cnt <- as.data.frame(t(counts(sub_dds, normalized=TRUE)))
  colnames(norm.cnt) <- paste0("norm_", names(markers_id))

  output <- output %>%
    bind_cols(cnt[output$sample_name, ]) %>% 
    bind_cols(norm.cnt[output$sample_name, ]) %>%
    arrange(dux4.rlogsum)

  write.xlsx(output, file=file_name,
             append=TRUE, row.names=FALSE,
             sheetName="dux4.group")
}

.makeFourBiomarkerHeatmap <- function(rlg, group,
                                      markers_id, file_name) {
  require(pheatmap)
  data <- assay(rlg[markers_id])
  rownames(data) <- names(markers_id)
  sample_name <- colnames(rlg)
  annotation_col <- data.frame(pheno_type=rlg$pheno_type,
                               DUX4_group=factor(group[sample_name]))
  rownames(annotation_col) <- sample_name                               
  pheatmap(data, treeheight_row=0,
           cellheight=16,
           annotation_col = annotation_col,
           file=file_name, slient=TRUE)                                          
}

.make53BiomarkerHeatmap <- function(rlg, group,
                                      markers_id, file_name) {
  require(pheatmap)
  data <- assay(rlg[markers_id])
  rownames(data) <- names(markers_id)
  sample_name <- colnames(rlg)
  annotation_col <- data.frame(pheno_type=rlg$pheno_type,
                               DUX4_group=factor(group[sample_name]))
  rownames(annotation_col) <- sample_name                               
  pheatmap(data, treeheight_row=0,
           cellheight=15,
           annotation_col = annotation_col,
           file=file_name, slient=TRUE)                                          
}

.makeFourBiomarkerScatterPlot <- function(dds, 
                              group, markers_id, 
                              call_threshold=2, file_name,
                              plot.tile_top = FALSE) {
  require(ggplot2)
  require(gridExtra)                               
  #' get local rlog, relaitve rlog, and make calls
  local_rlg <- .get_local_rlg(markers_id, dds)
  rownames(local_rlg) <- names(markers_id)
  centered <- rowMeans(local_rlg[, dds$pheno_type == "Control"])
  relative_rlg <- local_rlg - centered
  calls <- colSums(relative_rlg> call_threshold)  
  relative_logsum <- colSums(relative_rlg)
  melted <- as.data.frame(relative_rlg) %>% 
    rownames_to_column(var="gene_name") %>% 
    gather(key=sample_ID, value=rlogsum, -gene_name)
  
  #' prepate data for scatter plot
  sample_name <- colnames(dds)
  df <- tibble(
    sample_ID = sample_name,
    sum_rlog = relative_logsum[sample_name],
    DUX4_group = factor(group[sample_name]),
    calls = calls[sample_name]) %>%
    mutate(sample_ID = ifelse(grepl("b", sample_ID), paste0(sample_ID, "*"), sample_ID)) %>%
    arrange(sum_rlog) 
  sample_levels <- as.character(df$sample_ID)
  df$sample_ID <- factor(df$sample_ID, 
                         levels = sample_levels)

  #' scatter plot
  scatterPlot = ggplot(df, aes(x=sample_ID, y=sum_rlog,
                             color=DUX4_group)) +
    geom_point(aes(size=calls+1)) +
    #guides(size=FALSE) +
    scale_size(name = "number of calls", labels = c("0", "1", "2", "3", "4")) +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
          legend.justification=c(0,1), legend.position=c(0.02, 1),
          #legend.box.background = element_rect(colour = "black"),
          legend.text=element_text(size=9),
          legend.title=element_text(size=9)) +
    labs(y="sum of relative expression \n sum(rlog - average control)",
         x="sample name") 

  #' tileplot
  melted$sample_ID <- 
    factor(melted$sample_ID, levels=sample_levels)
  tilePlot <- ggplot(melted, aes(y=gene_name, x=sample_ID,
                                 fill=rlogsum)) +
    geom_tile() + 
    theme_minimal() +
    theme(legend.direction="horizontal",
          legend.position="top", 
          legend.text=element_text(size=8),
          legend.title=element_text(size=9),
          legend.justification=c(1, 0),
          axis.text.x=element_blank(), 
          axis.title.x=element_blank(),
          axis.title=element_blank(),  
          axis.text.y=element_text(size=8)) +
    scale_fill_gradient(low="white", high="steelblue", 
                        name="relative expression") +
    guides(fill = guide_colorbar(barwidth = 7, barheight = 0.7,
                                 title.position = "top", 
                                 title.hjust = 0.5))

  if (plot.tile_top) {                                 
    pdf(file_name, height=6, width=6)
    grid.arrange(tilePlot, scatterPlot, nrow=2, ncol=1, heights=c(1,3))
    dev.off()                
  }
  if (!plot.tile_top) {
    pdf(file_name, height=4.4, width=6.5)
    plot(scatterPlot)
    dev.off()
  }
  return(list(tilePlot=tilePlot, scatterPlot=scatterPlot))     

}


.getExpFromDDS <- function(gene_name, dds=year2.dds,
                           type=c("counts", "TPM", "RPKM")) {
  type <- match.arg(type)
  data <- tibble(gene_name=gene_name,
                 gene_id=get_ensembl(rse=dds, 
                                     inquiry_genename = gene_name)) 
  dds <- dds[data$gene_id]
  exp <- switch(type,
                counts = counts(dds, normalized=TRUE),
                TPM = assays(dds)[["TPM"]],
                RPKM = assays(dds)[["RPKM"]])
 
  data <- bind_cols(data, as.tibble(exp))
  data <- gather(data, sample_name, value, 
                 -gene_id, -gene_name)
  #' group = G1_C, G1_F, G2, and G3
  group <- substr(
    as.character(colData(dds)[data$sample_name, "pheno_type"]),
    start=1, stop=1)
  group <- paste0("G", colData(dds)[data$sample_name, "dux4.group"],
                  "_", group)
  data <- 
    add_column(data,
               DUX4_group=colData(dds)[data$sample_name, "dux4.group"],
               group=group)
}

.getExpFromDDS_by_cluster <- function(gene_name, dds=sanitized.dds,
                                     factor,
                           type=c("counts", "TPM", "RPKM")) {
  type <- match.arg(type)
  data <- tibble(gene_name=gene_name,
                 gene_id=get_ensembl(rse=dds, 
                                     inquiry_genename = gene_name)) 
  dds <- dds[data$gene_id]
  exp <- switch(type,
                counts = counts(dds, normalized=TRUE),
                TPM = assays(dds)[["TPM"]],
                RPKM = assays(dds)[["RPKM"]])
 
  data <- bind_cols(data, as.tibble(exp))
  data <- gather(data, sample_name, value, 
                 -gene_id, -gene_name)
  #' group = Control, Mild, Moderate, IG-High, Hot, Muscle-Low
  data <- 
    add_column(data,
               group=factor[data$sample_name])
}

.clean_covariate_for_corrgram <- function(time_var) {
  tidy_time <- time_var %>%
    select(DUX4.Score, ExtraCellular.Matrix.Score,
           Immune.Inflamm.Score, Immunoglobulin.Score,
           fat_rating, FF, STIR_rating,
           Pathlogy.Score, Fibrosis, Inflammation,
           Variability.in.Fiber.Size, 
           Extent.of.Central.Nucleation) %>%
    filter(!is.na(Inflammation)) %>%
    rename(fat_fraction=FF, Pathology.Score=Pathlogy.Score)

  #' turn factor to numerical
  tidy_time$Fibrosis <- 
    as.numeric(factor(tidy_time$Fibrosis,
      levels=c("Normal", "Mild", "Moderate", "Severe")))
  tidy_time$Variability.in.Fiber.Size <-
    as.numeric(factor(tidy_time$Variability.in.Fiber.Size,
      levels=c("Normal", "Mild", "Moderate", "Severe")))
  tidy_time$Extent.of.Central.Nucleation <-
    as.numeric(factor(tidy_time$Extent.of.Central.Nucleation,
       levels=c("Normal", "Mild", "Moderate", "Severe")))
  tidy_time
}

.adjust_time_2 <- function(time_2, relative_rlg) {
  #' (1) STIR_rating: round up
  adjust_time_2 <- time_2 %>%
    mutate(STIR_rating=round(STIR_rating))
  #' (2) fat_rating: 32-0015 -> 0
  adjust_time_2[adjust_time_2$sample_name == "32-0015b", "fat_rating"] <- 0  

  #' (3) recap and rename: RNA.score -> dux2.score, ii.score...
  sample_name <- as.character(adjust_time_2$sample_name)
  markers_score <- lapply(relative_rlg, function(x) x - min(x))
  markers_score <- as.data.frame(do.call(cbind, markers_score))
  rownames(markers_score)[rownames(markers_score)=="32-0002b1"] <- "32-0002b"

  adjust_time_2 <- adjust_time_2 %>%
    add_column(
      DUX4.Score              = markers_score[sample_name, "dux4"],
      ExtraCellular.Matrix.Score    = markers_score[sample_name, "extracellular_matrix"],
      Immune.Inflamm.Score       = markers_score[sample_name, "inflamm"],
      Immunoglobulin.Score  = markers_score[sample_name, "immunoglobulin"]) 
  
  return(adjust_time_2)
}

panel.corrgram.2 <-
    function(x, y, z, subscripts, at=pretty(z), scale=0.8, ...) {
        require("grid", quietly=TRUE)
        x <- as.numeric(x)[subscripts]
        y <- as.numeric(y)[subscripts]
        z <- as.numeric(z)[subscripts]
        zcal <- level.colors(z, at=at, ...)
        for (i in seq(along=z)) {
            lims <- range(0, z[i])
            tval <- 2 * base::pi *
                seq(from=lims[1], to=lims[2], by=0.01)
            grid.polygon(x=x[i] + 0.5 * scale * c(0, sin(tval)),
                         y=y[i] + 0.5 * scale * c(0, cos(tval)),
                         default.units="native",
                         gp=gpar(fill=zcal[i]))
            grid.circle(x=x[i], y=y[i], r=0.5*scale,
                        default.units="native")
        }
    }

#'
#' Tools for Part 2 analysis
#'
.do_goseq <- function(universal, selected_genes, p_value=0.01, 
                      return.DEInCat=FALSE, dds=NULL) {
    require(goseq)
    library(org.Hs.eg.db)
    library(GO.db)
    library(hg38.HomoSapiens.Gencode.v24)

    txsByGene <- transcriptsBy(hg38.HomoSapiens.Gencode.v24, by="gene")
    names(txsByGene) <- 
        sapply(strsplit(names(txsByGene), ".", fixed=TRUE), "[[", 1)
    lengthData <- median(width(txsByGene))
    
    isDEGs <- as.integer(universal %in% selected_genes)
    names(isDEGs) <- universal
    bias.data <- lengthData[names(isDEGs)]
    pwf <- nullp(isDEGs, bias.data=bias.data, plot.fit=FALSE)
    GO.BP <- goseq(pwf, "hg38", "ensGene", test.cats=c("GO:BP"))
    
    enriched.BP <- GO.BP %>%
      mutate(padj = p.adjust(over_represented_pvalue, method="BH")) %>%
      filter(padj < p_value)
    
    if (return.DEInCat & !is.null(dds)) {
      cat_genes <- lapply(enriched.BP$category, function(GOID) {
        cat_genes <- .cat2DEgenes(GOID, pwf=pwf)
        cat_genename <- .mapID2GeneName(dds, 
                                        id=cat_genes, clean=TRUE)
        paste(cat_genename, collapse=",")
      })  
      enriched.BP <- add_column(enriched.BP, DEInCat=cat_genes)
    }

    return(enriched.BP)
}

.cat2DEgenes <- function(GOID, pwf) {
    gene2cat <- getgo(rownames(pwf), "hg38", "ensGene", fetch.cats = "GO:BP")
    names(gene2cat) <- rownames(pwf)
    cat2gene <- goseq:::reversemapping(gene2cat)
    #' sanity check
    doesIDexist <- GOID %in% names(cat2gene)
    if (!all(doesIDexist)) stop("GOID is not found")
    sig_gene <- rownames(pwf)[pwf$DEgenes==1]
    geneInGOID <- cat2gene[[GOID]]
    sig_gene[sig_gene %in% geneInGOID]
}

.mapID2GeneName <- function(dds, id, clean=TRUE) {
    if (clean)  {
      rownames(dds) <- sapply(strsplit(rownames(dds), ".", fixed=TRUE),
                    "[[", 1)
    }
    rowData(dds[id])$gene_name
}

#'
#' Part 2 analysis; might support some part 1 analysis
#'

#' calculate local rlog
.get_local_rlg <- function(markers_id, dds) {
  flag <- markers_id %in% rownames(dds)
  if (!all(flag)) {
    stop(markers_id[!flag], " are not on the gene list of dds.")
  }
  markers_dds <- dds[markers_id]
  markers_rlg <- rlog(markers_dds, fitType="mean", blind=TRUE) # replace local
  local_rlg <- assay(markers_rlg)
}


#' calculate sum of relative local rlog
.get_sum_relative_local_rlg <- function(markers_id, dds) {
  local_rlg <- .get_local_rlg(markers_id, dds)
  centered <- rowMeans(local_rlg[, dds$pheno_type == "Control"])
  relative_rlg <- local_rlg - centered
  colSums(relative_rlg) 
}


#' get zscore of rlog expression
