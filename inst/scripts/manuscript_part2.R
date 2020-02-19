#'
#' load library
#'
library(DESeq2)
library(xlsx)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(gridExtra)
library(BiocParallel)
multi_param <- MulticoreParam(worker=4)
register(multi_param, default=TRUE)

pkg_dir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsy.all"
fig_dir <- file.path(pkg_dir, "manuscript", "figures")
table_dir <- file.path(pkg_dir, "manuscript", "tables")
sup_dir <- file.path(pkg_dir, "manuscript", "sup_figs")


#'
#' load common data: time_1, time_1, sanitized.dds, sanitized.rlg
#'                   year2.rlg, year2.dds
#'
source(file.path(pkg_dir, "scripts", "tools.R"))
source(file.path(pkg_dir, "scripts", "manuscript_tools.R"))

load(file.path(pkg_dir, "public_data", "sanitized.dds.rda"))
load(file.path(pkg_dir, "public_data", "sanitized.rlg.rda"))

load(file.path(pkg_dir, "public_data", "year2.dds.rda")) # from manuscript_part1.R
load(file.path(pkg_dir, "public_data", "year2.rlg.rda")) # from manuscript_part1.R

#'
#' Normalized visit I and II samples all together and get local rlog
#'
markers_list <- list(dux4 = c("LEUTX", "KHDC1L", "PRAMEF2", "TRIM43"),
                     extracellular_matrix = c("PLA2G2A", "COL19A1", "COMP", "COL1A1"), # remove SFPR2
                     inflamm = c("CCL18", "CCL13" ,"C6", "C7"), # remove CCL19
                     cell_cycle = c("CCNA1", "CDKN1A", "CDKN2A"), # newly added
                     immunoglobulin = c("IGHA1", "IGHG4", "IGHGP")) #delete IGKV3-11                    
markers_id_list <- lapply(markers_list, get_ensembl, rse=sanitized.dds)
relative_rlg <- lapply(markers_id_list, .get_sum_relative_local_rlg,
                       dds = sanitized.dds) 
scores <- as.data.frame(do.call(cbind, relative_rlg)) %>%
  rename(dux4.rlogsum=dux4, ecm.rlogsum=extracellular_matrix, inflamm.rlogsum=inflamm,
         cellcycle.rlogsum=cell_cycle,
         img.rlogsum=immunoglobulin) %>%
  rownames_to_column(var="sample_name") %>%
  mutate(pheno_type = sanitized.dds[, sample_name]$pheno_type,
         paired = sanitized.dds[, sample_name]$paired,
         patient_id = sanitized.dds[, sample_name]$patient_id,
         visit = sanitized.dds[, sample_name]$visit) %>%
  mutate(visit = factor(ifelse(visit=="I", "initial", "follow-up"), levels=c("initial", "follow-up"))) %>%
  dplyr::filter(pheno_type == "FSHD") %>%
  dplyr::filter(!is.na(paired)) %>%
  dplyr::filter(patient_id != "32-0008") 

#'######################################
#' (1) barplot for changes over two years
#'#######################################
#' DUX4 marker score
tmp <- scores %>% dplyr::filter(visit=="initial") %>% arrange(dux4.rlogsum)
patient_id <- as.character(tmp$patient_id)
scores$patient_id <- factor(scores$patient_id, levels=patient_id)
gg_dux4 <- ggplot(scores, aes(y=dux4.rlogsum, x=patient_id, fill=visit)) +
  geom_bar(stat="identity", position=position_dodge(0.7)) + 
  theme_minimal() + 
  labs(y="DUX4 score") +
  scale_fill_manual(values=c("grey75", "grey25")) +
  #scale_fill_brewer(palette="Greys") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=9),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=8),
        legend.justification=c(0,1), legend.position=c(0, 1),
        legend.key.size = unit(0.8, "line"))
pdf(file.path(fig_dir, "comb_dux4_score_barplot.pdf"), width=4, height=2.5)
plot(gg_dux4)
dev.off()

#' DUX4 scattor plot and correlation test for initial and follow-up visits
tmp <- scores %>% dplyr::select(patient_id, visit, dux4.rlogsum) %>%
  tidyr::spread(key=visit, -patient_id) %>% 
  dplyr::rename(follow_up = "follow-up")

#' ICC and correlation
irr2 <- irr::icc(tmp[, c(2,3)], type="consistency", unit="average", model="twoway") # 0.626
dux4_cor <- cor(tmp$initial, tmp$follow_up) # pearson

gg_dux4_scatter <- ggplot(tmp, aes(x=initial, y=follow_up)) +
  geom_point() +
  geom_text(x=0, y=15, hjust=0, 
            size=3, label=paste0("correlation = ", format(dux4_cor, digit=3))) +
  geom_smooth(method="lm") +
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5, size=12)) +
  labs(title="DUX4 score: initial vs. follow-up", y="follow-up visit", x="initial visit")
pdf(file.path(fig_dir, "comb_dux4_score_scatterplot.pdf"), width=4, height=2.5)
plot(gg_dux4_scatter)
dev.off()

#' ecm
tmp <- scores %>% dplyr::filter(visit=="initial") %>% arrange(ecm.rlogsum)
patient_id <- as.character(tmp$patient_id)
scores$patient_id <- factor(scores$patient_id, levels=patient_id)
gg_ecm <- ggplot(scores, aes(y=ecm.rlogsum, x=patient_id, fill=visit)) +
  geom_bar(stat="identity", position=position_dodge(0.7)) + 
  theme_minimal() + 
  labs(y="Extracellular matrix score") +
  scale_fill_manual(values=c("#E69F00", "#999999")) +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_text(size=9),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=8),
        legend.justification=c(0,1), legend.position=c(0, 1),
        legend.key.size = unit(0.8, "line"))

#' infl
tmp <- scores %>% dplyr::filter(visit=="initial") %>% arrange(inflamm.rlogsum)
patient_id <- as.character(tmp$patient_id)
scores$patient_id <- factor(scores$patient_id, levels=patient_id)
gg_inflam <- ggplot(scores, aes(y=inflamm.rlogsum, x=patient_id, fill=visit)) +
  geom_bar(stat="identity", position=position_dodge(0.7)) + 
  theme_minimal() + 
  labs(y="Inflammatory response score") +
  scale_fill_manual(values=c("#E69F00", "#999999")) +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_text(size=9),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=8),
        legend.justification=c(0,1), legend.position="none",
        legend.key.size = unit(0.8, "line"))

#' Ig
tmp <- scores %>% dplyr::filter(visit=="initial") %>% arrange(img.rlogsum)
patient_id <- as.character(tmp$patient_id)
scores$patient_id <- factor(scores$patient_id, levels=patient_id)
gg_img <- ggplot(scores, aes(y=img.rlogsum, x=patient_id, fill=visit)) +
  geom_bar(stat="identity", position=position_dodge(0.7)) + 
  theme_minimal() + 
  labs(y="Immunoglobulins score") +
  scale_fill_manual(values=c("#E69F00", "#999999")) +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_text(size=9),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=8),
        legend.justification=c(0,1), legend.position="none",
        legend.key.size = unit(0.8, "line"))

#' Cell cycle
tmp <- scores %>% dplyr::filter(visit=="initial") %>% arrange(cellcycle.rlogsum)
patient_id <- as.character(tmp$patient_id)
scores$patient_id <- factor(scores$patient_id, levels=patient_id)
gg_cycle <- ggplot(scores, aes(y=cellcycle.rlogsum, x=patient_id, fill=visit)) +
  geom_bar(stat="identity", position=position_dodge(0.7)) + 
  theme_minimal() + 
  labs(y="Cell cycle score") +
  scale_fill_manual(values=c("#E69F00", "#999999")) +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_text(size=9),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=8),
        legend.justification=c(0,1), legend.position="none",
        legend.key.size = unit(0.8, "line"))

library(gridExtra)
pdf(file.path(fig_dir, "comb_biomarker_scores_Barplot.pdf"), width=8, height=6)
grid.arrange(gg_ecm, gg_inflam, gg_img, gg_cycle, nrow=2, ncol=2)
dev.off()

#'##########################################
#' (2) delta score (follow-up - initial visit)
############################################
load(file.path(pkg_dir, "public_data", "mri_pathology.rda"))

year1_scores <- scores %>% dplyr::filter(visit == "initial")
year2_scores <- scores %>% dplyr::filter(visit == "follow-up")
tmp <- inner_join(year1_scores, year2_scores, by="patient_id",
                  suffix=c(".I", ".II"))
delta_scores <- data.frame(patient_id = as.character(tmp$patient_id),
                           dux4 = tmp$dux4.rlogsum.II - tmp$dux4.rlogsum.I,
                           ecm = tmp$ecm.rlogsum.II - tmp$ecm.rlogsum.I,
                           inflamm = tmp$inflamm.rlogsum.II - tmp$inflamm.rlogsum.I,
                           img = tmp$img.rlogsum.II - tmp$img.rlogsum.I,
                           cycle = tmp$cellcycle.rlogsum.II - tmp$cellcycle.rlogsum.I)   

comb_time <- left_join(mri_pathology$time_1, mri_pathology$time_2, by="patient_id",
                       suffix=c(".I", ".II")) %>%
  dplyr::filter(!is.na(paired.I)) 

delta_clinic <- comb_time %>%
  mutate(pathology = Pathology.Score.II - Pathology.Score.I) %>%
  mutate(STIR = STIR_rating.II - STIR_rating.II,
         T1 = T1_rating.II - T1_rating.I)    %>%
  select(patient_id, pathology, STIR, T1) 

comb_delta <- full_join(delta_scores, delta_clinic, by="patient_id") %>%
  dplyr::filter(!is.na(pathology)) %>%
  dplyr::filter(!is.na(dux4)) %>%
  dplyr::select(patient_id, dux4, pathology) %>%
  dplyr::rename(DUX4.Score=dux4) %>%
  arrange(DUX4.Score) %>%
  mutate(patient_id = as.character(patient_id)) %>%
  mutate(patient_id = factor(patient_id, levels=as.character(patient_id)))

melt_scores <- gather(comb_delta, "features", "delta_scores", -patient_id)
gg <- ggplot(melt_scores, aes(y=delta_scores, x=patient_id, fill=features,
                              group=features)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.7) +
  theme_minimal() +
  labs(y=expression(paste(Delta, "score: follow-up - initial")) , 
       title="Difference between two visits: pathology and DUX4 scores" ) +
  theme(axis.title.x=element_blank(), plot.title=element_text(hjust=0.5, size=11),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.key.size = unit(0.8, "line"),
        legend.justification=c(0,1), legend.position="bottom") +
  #scale_fill_brewer(palette = "Greys") 
  scale_fill_manual(values=c("grey30", "grey80"))
  #geom_rect(aes(xmin=7.5, xmax=8.5, ymin=-15, ymax=15), fill="grey", alpha=0)    
  
pdf(file.path(sup_dir, "comb_path_dux4_delta_score.pdf"), width=8, height=3)
plot(gg)
dev.off()  

########################################################
#'
#' (3) Pathology/MRI score between year 1 and year 2
#' take out 32-0008 and (32-0008b and 32-0016b and 32-0010 don't have histology)
########################################################
path_score <- bind_rows(mri_pathology$time_1, mri_pathology$time_2) %>%
  mutate(visit = ifelse(visit=="I", "initial", "follow-up")) %>%
  mutate(visit = factor(visit, levels=c("initial", "follow-up"))) %>%
  dplyr::filter(!is.na(paired)) %>% 
  dplyr::filter(!is.na(Pathology.Score)) %>%
  dplyr::filter(!patient_id %in% c("32-0008", "32-0010", "32-0016")) 

tmp <- path_score %>% dplyr::filter(visit=="initial") %>%
  arrange(Pathology.Score) %>% 
  mutate(patient_id=as.character(patient_id)) %>%
  dplyr::pull(patient_id)

path_score <- path_score %>%
  mutate(patient_id = factor(patient_id, levels=tmp))

#' median of pathology
path_score %>% group_by(visit) %>%
  summarize(median(Pathology.Score))

gg0 <- ggplot(path_score, aes(x=patient_id, y=Pathology.Score, 
                           group=visit, fill=visit)) +
    geom_bar(stat="identity", position=position_dodge(0.7)) + 
  theme_minimal() + 
  labs(y="Pathology score") +
  scale_fill_manual(values=c("#E69F00", "#999999")) +
  theme(axis.title.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=8),
        legend.justification=c(0,1), legend.position=c(0, 1),
        legend.key.size = unit(0.8, "line"))
pdf(file.path(fig_dir, "comb_path_Barplot.pdf"), width=4, height=2.5)
plot(gg0)
dev.off()                

#' inflammation
gg1 <- ggplot(path_score, aes(x=patient_id, y=Inflammation, 
                           group=visit, fill=visit)) +
    geom_bar(stat="identity", position=position_dodge(0.7)) + 
  theme_minimal() + 
  labs(y="Inflammation") +
  scale_fill_manual(values=c("#E69F00", "#999999")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.justification=c(0,1), legend.position="none")   
#' fiber size
gg2 <- ggplot(path_score, aes(x=patient_id, y=Variability.in.Fiber.Size, 
                           group=visit, fill=visit)) +
    geom_bar(stat="identity", position=position_dodge(0.7)) + 
  theme_minimal() + 
  labs(y="Variability in fiber size") +
  scale_fill_manual(values=c("#E69F00", "#999999")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.justification=c(0,1), legend.position="none")   
#' extent of central nucleation
gg3 <- ggplot(path_score, aes(x=patient_id, y=Extent.of.Central.Nucleation, 
                           group=visit, fill=visit)) +
    geom_bar(stat="identity", position=position_dodge(0.7)) + 
  theme_minimal() + 
  labs(y="Extent of central nucleation") +
  scale_fill_manual(values=c("#E69F00", "#999999")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.justification=c(0,1), legend.position="none")   
#' necrosis regeneration and inflammation
gg4 <- ggplot(path_score, aes(x=patient_id, y=Necrosis.Regeneration.Inflammation, 
                           group=visit, fill=visit)) +
    geom_bar(stat="identity", position=position_dodge(0.7)) + 
  theme_minimal() + 
  labs(y="Necrosis regeneration \ inflamm") +
  scale_fill_manual(values=c("#E69F00", "#999999")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.justification=c(0,1), legend.position="none")   
#' fibrosis
gg5 <- ggplot(path_score, aes(x=patient_id, y=Fibrosis, 
                           group=visit, fill=visit)) +
    geom_bar(stat="identity", position=position_dodge(0.7)) + 
  theme_minimal() + 
  labs(y="Fibrosis") +
  scale_fill_manual(values=c("#E69F00", "#999999")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.justification=c(0,1), legend.position="none")  

pdf(file.path(fig_dir, "comb_patholoby_scores_barplot.pdf"), width=10, height=8)
grid.arrange(gg0, gg1, gg2, gg3, gg4, gg5, nrow=3, ncol=2)
dev.off()
                                                      

# density plot of pathology score
path_sum <- path_score %>% group_by(visit) %>% 
  summarize(median=median(Pathology.Score),
            mean=mean(Pathology.Score),
            sd=sd(Pathology.Score)) 

gg <- ggplot(path_score, aes(x=Pathology.Score, color=visit, fill=visit)) +
    geom_density(alpha=0.5) +
    #geom_area(aes(y=..density..), alpha=0.5, stat="bin") +
    scale_fill_manual(values=c("grey75", "grey25")) +
    scale_color_manual(values=c("grey75", "grey25")) +
    theme_bw() +
    theme(legend.justification=c(1,1), legend.position=c(0.99, 0.99), 
          legend.key.size = unit(0.8, "line"),
          #legend.title=element_text(size=9),
          legend.title=element_blank(),
          plot.title = element_text(hjust=0.5, size=10),
          axis.title = element_text(size=9),
          axis.text=element_text(size=7)) +
    labs(title="Density of pathology score", x="Pathloby score")      
gg <- gg + geom_vline(data=path_sum, aes(xintercept=median, color=visit),
             linetype="dashed", show.legend=FALSE) 

pdf(file.path(fig_dir, "comb_pathology_density.pdf"), width=3.5, height=2.5)
plot(gg)
dev.off()        

#'##################################################
#' (4) DUX4 vs STIR+/All FSHD/control, other markers vs. STIR+/FSHD/Control
#' append markers scores to path_score (no need to exclude non-paried samples)
#'#################################################
path_score <- bind_rows(mri_pathology$time_1, mri_pathology$time_2) %>%
  dplyr::filter(!is.na(STIR_rating)) %>%
  mutate(STIR_status = ifelse(STIR_rating > 0, "STIR+", "MRI normal")) %>%
  mutate(STIR_status = ifelse((STIR_rating == 0 & T1_rating > 0), NA, STIR_status)) %>%
  select(-visit, -patient_id, -paired, -Event.Name)

scores <- as.data.frame(do.call(cbind, relative_rlg)) %>%
  rename(dux4.rlogsum=dux4, ecm.rlogsum=extracellular_matrix, inflamm.rlogsum=inflamm,
         cellcycle.rlogsum=cell_cycle,
         img.rlogsum=immunoglobulin) %>%
  rownames_to_column(var="sample_name") %>%
  mutate(pheno_type = sanitized.dds[, sample_name]$pheno_type,
         paired = sanitized.dds[, sample_name]$paired,
         patient_id = sanitized.dds[, sample_name]$patient_id,
         visit = sanitized.dds[, sample_name]$visit) %>%
  mutate(visit = factor(ifelse(visit=="I", "initial", "follow-up"), levels=c("initial", "follow-up")))

path_score$sample_name[path_score$sample == "32-0002b"] <- "32-0002b1" # match with pathology sample id  
marker_path_scores <- left_join(scores, path_score, by="sample_name") %>%
  select(-paired)

save(marker_path_scores, file=file.path(pkg_dir, "public_data", "marker_path_scores.rda"))


#marker_score <- bind_rows(year1_scores, year2_scores) %>%
  #dplyr::filter(!is.na(paired)) %>%
#  dplyr::filter(patient_id != "32-0008") 

#new_time <- time_var %>% 
#  dplyr::filter(!is.na(paired)) %>%
#  dplyr::filter(patient_id != "32-0008") %>%
#  dplyr::filter(!is.na(Inflammation)) %>%
#  filter(!(STIR_rating == 0 & fat_rating > 0)) %>% # remove T1>0 and T2=0 (two0)
#  mutate(STIR_status = ifelse(STIR_rating > 0, "STIR+", "MRI Normal"))  %>% 
#  mutate(sample_name = as.character(sample_name))
#new_time$sample_name[new_time$sample_name == "32-0002b"] <- "32-0002b1"

#' append STIR status to marker_score
#sum(new_time$sample_name %in% marker_score$sample_name)
#marker_score <-  marker_score %>%
#  left_join(new_time, by="sample_name")

.mutate_row <- function(self) {
  stir <- self %>% dplyr::filter(STIR_status == "STIR+") %>%
    mutate(status = "STIR+")
  bind_rows(self, stir)  
}

tmp <- marker_path_scores %>% 
  mutate(status=ifelse(pheno_type=="Control", "Control", "Total")) %>%
  group_by(visit) %>% .mutate_row(.) %>% ungroup(.) %>%
  mutate(group=paste0(visit, "_", status)) %>%
  mutate(group=factor(group, levels=c("initial_Control", "initial_Total", "initial_STIR+",
                                      "follow-up_Total", "follow-up_STIR+"))) %>%
  mutate(visit = factor(visit, levels=c("initial", "follow-up")))                                      

gg <- ggplot(tmp, aes(x=group, y=dux4.rlogsum, group=group)) +
  geom_boxplot(aes(color=visit), width=0.5) +
  #geom_jitter(width=0.2, size=1) +
  #scale_fill_manual(values=c("#E69F00", "#999999")) +
  #scale_color_manual(values=c("#E69F00", "#999999")) +
  scale_fill_manual(values=c("grey75", "grey25")) +
  scale_color_manual(values=c("grey75", "grey25")) +
  theme_minimal() +
  geom_vline(xintercept=3.5, linetype="dashed", color = "grey") +
  theme(#legend.position=c(0.01,0.99),
        legend.position="right",
        legend.justification=c(0,1), 
        #legend.key.size = unit(0.7, "line"),
        plot.title = element_text(hjust = 0.5, size=10),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=8),
        axis.title.y=element_text(size=10)) + 
  labs(y="Marker score", title="DUX4") +      
  scale_x_discrete(labels=c("Control", "All\nFSHD", "STIR+\nFSHD", "All\nFSHD", "STIR+\nFSHD")) +
  coord_cartesian(ylim=c(-0.05, 23))

pdf(file.path(fig_dir, "comb_dux4_T2.pdf"), width=3.6, height=2.4)
plot(gg)
dev.off()

tmp2 <- tmp %>% 
  select(sample_name, ecm.rlogsum, inflamm.rlogsum, img.rlogsum, cellcycle.rlogsum, group, visit) %>%
  gather(markers, scores, -sample_name, -group, -visit) %>%
  mutate(markers = factor(markers, levels=c("ecm.rlogsum", "inflamm.rlogsum",
                                            "img.rlogsum", "cellcycle.rlogsum")))

labels <- c(ecm.rlogsum="Extracellular matrix", 
            inflamm.rlogsum="Inflamm response", 
            img.rlogsum="Immunoglobulin",
            cellcycle.rlogsum="Cell cycle")

gg <- ggplot(tmp2, aes(x=group, y=scores, group=group)) + 
    geom_boxplot(aes(color=visit), width=0.5) +
  #scale_fill_manual(values=c("#E69F00", "#999999")) +
  #scale_color_manual(values=c("#E69F00", "#999999")) +
  scale_fill_manual(values=c("grey75", "grey25")) +
  scale_color_manual(values=c("grey75", "grey25")) +
  theme_minimal() +
  labs(y="Marker score") +
  geom_vline(xintercept=3.5, linetype="dashed", color = "grey") +
  facet_wrap( ~ markers, nrow=2, scale="free", labeller=labeller(markers=labels)) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.x = element_text(size=8)) +
  scale_x_discrete(labels=c("Control", "All\nFSHD", "STIR+\nFSHD", "All\nFSHD", "STIR+\nFSHD"))   

pdf(file.path(fig_dir, "comb_markers_T2.pdf"), width=4.5, height=4)
plot(gg)
dev.off()


#'
#' fixing Suppl Table 1:
#'
load(file.path(pkg_dir, "public_data", "marker_path_scores.rda")
tmp <- marker_path_scores %>% 
  dplyr::filter(pheno_type=="FSHD") %>%
  dplyr::select(patient_id, visit, dux4.rlogsum) %>%
  spread(key=visit, value=dux4.rlogsum)

library(xlsx)
clinic <- read.xlsx(file.path(table_dir, "suppl_table_1_mri_pathology_dux4_scores.xlsx"), sheetIndex=1) %>%
  dplyr::filter(!is.na(Subject..), !is.na(Age)) %>%
  rename(patient_id = "Subject..") %>% 
  dplyr::select(patient_id, Age, Baseline.Muscle.Pathology) %>%
  left_join(tmp, by="patient_id")

year2_dux4_group <- as.data.frame(get(load(file.path(pkg_dir, "public_data", "year2_dux4_score.rda")))) %>%
  dplyr::select(sample_name, dux4.group) %>%
  mutate(patient_id = str_replace(as.character(sample_name), "b", "")) %>%
  dplyr::select(-sample_name) %>%
  rename(dux4.rank.followup=dux4.group) 

clinic <- clinic %>% left_join(year2_dux4_group, by="patient_id") %>%
  mutate(initial=format(initial, digit=4), `follow-up`=format(`follow-up`, digit=4))
  
write.csv(clinic, file=file.path(table_dir, "mri_pathology_dux4_scores.csv"))  

#'
#' Fixing Suppl Table 2 
#'
path <- read.xlsx(file.path(table_dir, "supple_table_2_pathology_feature_scores.xlsx"), 
                  sheetIndex=1, stringsAsFactors=FALSE) %>%
  mutate(new_id = get_alt_name(Record.ID) )                  
write.csv(path, file=file.path(table_dir, "patholog_feature_scores_newname.csv"))
