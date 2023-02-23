##Load the necessary packages
library(rlang)
library(tidyverse)
library(BiocManager)
library(clusterProfiler)
library(GGally)
library(ggrepel)
library(limma)
library(RColorBrewer)
library(msigdbr)
library(GOSemSim)
library(enrichplot)
library(ggalt)
library(ggExtra)
library(ggpubr)
library(cowplot)
library(DOSE)
library(rje)
library(DEqMS)
library(BiocParallel)
library(ggthemes)
library(Hmisc)
library(PTXQC)
library(Biostrings)
library(matrixStats)
library(data.table)
library(readxl)
library(org.Hs.eg.db)
library(ggh4x)
library(iq)

#### set working directory !!!!!!!
setwd("~ path to your working directory")
getwd()


### raw data and DIA-NN output tables are available via PRIDE Id PXD038915 !

# import protein groups table
precursors <- read_tsv(file = "report.pr_matrix.tsv") # intensities with no channel information
colnames(precursors) <- c("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description","Proteotypic","Stripped.Sequence","Modified.Sequence","Precursor.Charge","Precursor.Id",
                              paste0("2h_R",seq(1:3)), paste0("4h_R",seq(1:3)), paste0("6h_R",seq(1:3)), paste0("9h_R",seq(1:3)), paste0("24h_R",seq(1:3)))

precursors <- precursors[!grepl("Cont_", precursors$Protein.Ids),] # remove contaminants

pr_l <- precursors %>% gather(contains(c("_R")), key = "Experiment", value = "LFQ")

pr_l$sample = substr(pr_l$Experiment,1,nchar(pr_l$Experiment)-3)

ident_prec <- pr_l[!is.na(pr_l$LFQ),] # remove missing values from long protein groups df
ident_prec$n <- rep(1)
ident_prec_n <- aggregate(data = ident_prec, `n` ~ `Experiment` + `Protein.Ids` + `Genes` + `sample`, FUN = sum, na.rm = TRUE)
n_iprec <- aggregate(data = ident_prec_n, `n` ~ Experiment + `sample`, FUN = sum, na.rm = TRUE)

# plot numberof LFQ quantified proteins
iprec_count_mean <- aggregate(data = n_iprec, `n` ~ sample, FUN = mean, na.rm = TRUE)

# plot number of identified precursors (excluding different SILAC labels)
ident_prec <- ggplot() + 
  geom_col(iprec_count_mean, mapping = aes(x = as.factor(`sample`), fill = `sample`, y = `n`), position = "dodge", show.legend = FALSE, alpha = 0.7) +
  geom_jitter(n_iprec, mapping = aes(x = as.factor(`sample`), group = `sample`, y = `n`), width = 0.25, show.legend = FALSE, size = 3, alpha = 0.6) +
  geom_text(iprec_count_mean, mapping = aes(x = as.factor(`sample`), group = `sample`, label = round(`n`, digits = 0), y = `n`), vjust = 0.5, size = 5, 
            color = "black", position = position_dodge(width = 0.95), angle = 90, hjust = 1.5) + 
  ylab("identified precursors") + xlab("") + theme_bw() + scale_x_discrete(limits = c("2h","4h","6h","9h","24h")) +
  theme(plot.title = element_text(color="black", size = 16, face= "bold", vjust = 0.5), axis.title.x = element_text(color="black", size=14, face= "bold"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color="black", size = 14), axis.text.y = element_text(color="black", size=14), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(), axis.title.y = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 14, face = "bold"), legend.position = "right",
        legend.text = element_text(color="black", size=14), legend.title = element_text(color="black", size=14, face="bold"))

ident_prec

# import protein groups table
protein.groups <- read_tsv(file = "report.pg_matrix.tsv") # intensities with no channel information
colnames(protein.groups) <- c("Protein.Groups","Protein.Ids","Protein.Names","Genes","First.Protein.Description",
                             paste0("2h_R",seq(1:3)), paste0("4h_R",seq(1:3)), paste0("6h_R",seq(1:3)), paste0("9h_R",seq(1:3)), paste0("24h_R",seq(1:3)))

protein.groups <- protein.groups[!grepl("Cont_", protein.groups$Protein.Ids),] # remove contaminants

lfq_l <- protein.groups %>% gather(contains(c("_R")), key = "Experiment", value = "LFQ")

lfq_l$sample = substr(lfq_l$Experiment,1,nchar(lfq_l$Experiment)-3)

ident_prot <- lfq_l[!is.na(lfq_l$LFQ),] # remove missing values from long protein groups df
ident_prot$n <- rep(1)
ident_prot_n <- aggregate(data = ident_prot, `n` ~ `Experiment` + `Protein.Ids` + `Genes` + `sample`, FUN = sum, na.rm = TRUE)
n_iprot <- aggregate(data = ident_prot_n, `n` ~ Experiment + `sample`, FUN = sum, na.rm = TRUE)

# plot numberof LFQ quantified proteins
iprot_count_mean <- aggregate(data = n_iprot, `n` ~ sample, FUN = mean, na.rm = TRUE)

# plot number of identified protein groups (no SILAC ratio quantifiaction)
ident_prot <- ggplot() + 
  geom_col(iprot_count_mean, mapping = aes(x = as.factor(`sample`), fill = `sample`, y = `n`), position = "dodge", show.legend = FALSE, alpha = 0.7) +
  geom_jitter(n_iprot, mapping = aes(x = as.factor(`sample`), group = `sample`, y = `n`), width = 0.25, show.legend = FALSE, size = 3, alpha = 0.6) +
  geom_text(iprot_count_mean, mapping = aes(x = as.factor(`sample`), group = `sample`, label = round(`n`, digits = 0), y = `n`), vjust = 0.5, size = 5, 
            color = "black", position = position_dodge(width = 0.95), angle = 90, hjust = 1.5) + #scale_x_discrete(breaks = c("2h","4h","6h","9h","24h")) +
  ylab("identified protein groups") + xlab("") + theme_bw() + scale_x_discrete(limits = c("2h","4h","6h","9h","24h")) +
  theme(plot.title = element_text(color="black", size = 16, face= "bold", vjust = 0.5), axis.title.x = element_text(color="black", size=14, face= "bold"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color="black", size = 14), axis.text.y = element_text(color="black", size=14), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(), axis.title.y = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 14, face = "bold"), legend.position = "right",
        legend.text = element_text(color="black", size=14), legend.title = element_text(color="black", size=14, face="bold"))

ident_prot

# channel specific MS2-based precursor translated quantity
p_translated <- read_tsv(file = "report.pr_matrix_channels_translated.tsv") # MS2-based quantification
p_translated <- p_translated[!grepl("Cont_", p_translated$Protein.Ids),] # remove contaminants

p_translated_L <- p_translated %>% dplyr::select(`Protein.Group`:`Precursor.Id`, contains("raw-L"))
colnames(p_translated_L) = c("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description","Proteotypic","Stripped.Sequence","Modified.Sequence","Precursor.Charge","Precursor.Id",
                             paste0("2h_L",seq(1:3)), paste0("4h_L",seq(1:3)), paste0("6h_L",seq(1:3)), paste0("9h_L",seq(1:3)), paste0("24h_L",seq(1:3)))

p_translated_M <- p_translated %>% dplyr::select(`Protein.Group`:`Precursor.Id`, contains("raw-M"))
colnames(p_translated_M) = c("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description","Proteotypic","Stripped.Sequence","Modified.Sequence","Precursor.Charge","Precursor.Id",
                             paste0("2h_M",seq(1:3)), paste0("4h_M",seq(1:3)), paste0("6h_M",seq(1:3)), paste0("9h_M",seq(1:3)), paste0("24h_M",seq(1:3)))

p_translated_H <- p_translated %>% dplyr::select(`Protein.Group`:`Precursor.Id`, contains("raw-H"))
colnames(p_translated_H) = c("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description","Proteotypic","Stripped.Sequence","Modified.Sequence","Precursor.Charge","Precursor.Id",
                             paste0("2h_H",seq(1:3)), paste0("4h_H",seq(1:3)), paste0("6h_H",seq(1:3)), paste0("9h_H",seq(1:3)), paste0("24h_H",seq(1:3)))

p_translated <- full_join(p_translated_L, p_translated_M) %>% full_join(p_translated_H)

pg_l <- p_translated %>% gather(contains(c("_L","_M","_H")), key = "Experiment", value = "precursor_translated")
pg_l$precursor_translated[pg_l$precursor_translated == 0] <- NA

SILAC_precursors <- pg_l[!is.na(pg_l$`precursor_translated`),] # remove 0 values for norm SILAC ratios

# add specified label info
prec_S <- SILAC_precursors %>% mutate(label = case_when(endsWith(Experiment, "_L1") ~ "L", endsWith(Experiment, "_L2") ~ "L", endsWith(Experiment, "_L3") ~ "L",
                                                        endsWith(Experiment, "_M1") ~ "M", endsWith(Experiment, "_M2") ~ "M", endsWith(Experiment, "_M3") ~ "M",
                                                        endsWith(Experiment, "_H1") ~ "H", endsWith(Experiment, "_H2") ~ "H", endsWith(Experiment, "_H3") ~ "H"))


# create table with number of distinct precursor numbers (excluding SILAC channels) for each protein group - used later for DEqMS
n_precursors <- prec_S[,c("Experiment","Precursor.Id","Protein.Ids","Genes")]
n_precursors$Experiment <- gsub("_L", "_R", n_precursors$Experiment)
n_precursors$Experiment <- gsub("_M", "_R", n_precursors$Experiment)
n_precursors$Experiment <- gsub("_H", "_R", n_precursors$Experiment)
n_precursors <- unique(n_precursors)
n_precursors$n <- rep(1)
prec_count <- aggregate(data = n_precursors, `n` ~ `Experiment` + `Protein.Ids` + `Genes`, FUN = sum, na.rm = TRUE)


write_tsv(prec_S[,c("Protein.Ids","Genes","Precursor.Id","Stripped.Sequence","Modified.Sequence","Experiment","precursor_translated")], "all_precursors_ldf_MS2.tsv") # save table of precurosrs for IQ package to use for LFQ calculation

process_long_format("all_precursors_ldf_MS2.tsv", output_filename = "iq_pg_proc_MS2.tsv", primary_id = "Protein.Ids", secondary_id = c("Precursor.Id","Stripped.Sequence","Modified.Sequence"),
                    annotation_col = c("Genes"), filter_double_less = NULL, intensity_col_sep = ";", sample_id = "Experiment", intensity_col = "precursor_translated") # iq package based  calculation of LFQ values for the separate SILAC channels

# import data
prot_lfq <- read_tsv("iq_pg_proc_MS2.tsv")

pg_l <- prot_lfq %>% gather(contains(c("_L","_M","_H")), key = "Experiment", value = "LFQ")

pg_l <- pg_l %>% mutate(label = case_when(endsWith(Experiment, "_L1") ~ "L", endsWith(Experiment, "_L2") ~ "L", endsWith(Experiment, "_L3") ~ "L",
                                          endsWith(Experiment, "_M1") ~ "M", endsWith(Experiment, "_M2") ~ "M", endsWith(Experiment, "_M3") ~ "M",
                                          endsWith(Experiment, "_H1") ~ "H", endsWith(Experiment, "_H2") ~ "H", endsWith(Experiment, "_H3") ~ "H"))

# sub channel info in experiment title back to replicate
pg_l$Experiment <- str_replace_all(pg_l$Experiment, "_L", "_R")
pg_l$Experiment <- str_replace_all(pg_l$Experiment, "_M", "_R")
pg_l$Experiment <- str_replace_all(pg_l$Experiment, "_H", "_R")

# split into SILAC channels
prot_L <- pg_l %>% filter(label == "L")
#prot_L <- pg_l %>% filter(stringr::str_detect(Experiment, '_L'))
colnames(prot_L) <- c("Protein.Ids","Genes","Experiment","int_L","label")
prot_M <- pg_l %>% filter(label == "M")
#prot_M <- pg_l %>% filter(stringr::str_detect(Experiment, '_M1|_M2|_M3'))
colnames(prot_M) <- c("Protein.Ids","Genes","Experiment","int_M","label")
prot_H <- pg_l %>% filter(label == "H")
#prot_H <- pg_l %>% filter(stringr::str_detect(Experiment, '_H'))
colnames(prot_H) <- c("Protein.Ids","Genes","Experiment","int_H","label")

prot_SF <- full_join(prot_L[,1:4], prot_M[,1:4]) %>% full_join(prot_H[,1:4])
prot_SF$`ratio_H/M` <- 2^(prot_SF$int_H) / 2^(prot_SF$int_M)
prot_SF$`ratio_H/L` <- 2^(prot_SF$int_H) / 2^(prot_SF$int_L)
prot_SF$`ratio_M/L` <- 2^(prot_SF$int_M) / 2^(prot_SF$int_L)

prot_SF$sample <- substr(prot_SF$Experiment,1,nchar(prot_SF$Experiment)-3)
prot_SF$perc <- substr(prot_SF$Experiment,1,nchar(prot_SF$Experiment)-3)


# aggregate SILAC ratios of precursors to median of proteins
prot_HM_agg <- aggregate(data = prot_SF, `ratio_H/M` ~ Experiment + Protein.Ids  + `Genes` + perc + sample, FUN = median, na.rm = TRUE)
prot_HM_agg$ratio <- rep("H/M")
colnames(prot_HM_agg) <- c("Experiment","Protein.Ids","Genes","perc","sample","num_ratio","ratio")

prot_HL_agg <- aggregate(data = prot_SF, `ratio_H/L` ~ Experiment + Protein.Ids  + `Genes` + perc + sample, FUN = median, na.rm = TRUE)
prot_HL_agg$ratio <- rep("H/L")
colnames(prot_HL_agg) <- c("Experiment","Protein.Ids","Genes","perc","sample","num_ratio","ratio")

prot_ML_agg <- aggregate(data = prot_SF, `ratio_M/L` ~ Experiment + Protein.Ids  + `Genes` + perc + sample, FUN = median, na.rm = TRUE)
prot_ML_agg$ratio <- rep("M/L")
colnames(prot_ML_agg) <- c("Experiment","Protein.Ids","Genes","perc","sample","num_ratio","ratio")

prot_HL_agg$ratio <- rep("H/L")
prot_HM_agg$ratio <- rep("H/M")
prot_ML_agg$ratio <- rep("M/L")

colnames(prot_HM_agg) <- c("Experiment","Protein.Ids","Genes","perc","sample","num_ratio","ratio")
colnames(prot_HL_agg) <- c("Experiment","Protein.Ids","Genes","perc","sample","num_ratio","ratio")
colnames(prot_ML_agg) <- c("Experiment","Protein.Ids","Genes","perc","sample","num_ratio","ratio")

# combine and arrange
prot_SILAC_agg <- rbind(prot_HM_agg, prot_HL_agg, prot_ML_agg)

# save and import
#write_tsv(prot_SILAC_agg, "aggregated_SILAC_ratios.tsv")
prot_SILAC_agg <- read_tsv(file = "aggregated_SILAC_ratios.tsv")
prot_SILAC_agg <- prot_SILAC_agg[!is.na(prot_SILAC_agg$num_ratio),] # remove missing values from long datble

prot_SILAC_agg$n <- rep(1)
prot_count <- aggregate(data = prot_SILAC_agg, `n` ~ `Experiment` + `perc` + `sample` + `ratio`, FUN = sum, na.rm = TRUE)

# select sample
prot_SILAC_agg_f <- prot_SILAC_agg %>% dplyr::filter(`ratio` == "H/M") # use H/M SILAC ratios for quant comparison of IFNg treated cells / BSA control cells
prot_count_f <- prot_count %>% dplyr::filter(`ratio` == "H/M") 

prot_count_mean <- aggregate(data = prot_count_f, `n` ~ perc + sample, FUN = mean, na.rm = TRUE)

# plot
quant_prot <- ggplot() + 
  geom_col(prot_count_mean, mapping = aes(x = `sample`, fill = `sample`, y = `n`), position = "dodge", show.legend = FALSE, alpha = 0.7) +
  geom_jitter(prot_count_f, mapping = aes(x = `sample`, group = `sample`, y = `n`), width = 0.25, show.legend = FALSE, size = 3, alpha = 0.6) +
  geom_text(prot_count_mean, mapping = aes(x = `sample`, group = `sample`, label = round(`n`, digits = 0), y = `n`), vjust = 0.5, size = 5, 
            color = "black", position = position_dodge(width = 0.95), angle = 90, hjust = 1.5) + scale_x_discrete(limits = c("2h","4h","6h","9h","24h")) +
  ylab("quant. protein groups") + xlab("IFNg treatment") + theme_bw() + #ylim(0,7000) + #ggtitle("quantified protein groups") + #scale_fill_manual(values = c("dodgerblue3","firebrick3")) +
  theme(plot.title = element_text(color="black", size = 16, face= "bold", vjust = 0.5), axis.title.x = element_text(color="black", size=14, face= "bold"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color="black", size = 14), axis.text.y = element_text(color="black", size=14), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(), axis.title.y = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 14, face = "bold"), legend.position = "right",
        legend.text = element_text(color="black", size=14), legend.title = element_text(color="black", size=14, face="bold"))

quant_prot

# median norm
prot_silac_median <- aggregate(data = prot_SILAC_agg_f, `num_ratio` ~ perc + sample + ratio + Experiment, FUN = median, na.rm = TRUE) # median normalization of samples by subtraction of median
prot_silac_median <- rename(prot_silac_median, c('num_ratio' = 'median'))
prot_SILAC_agg_f <- full_join(prot_SILAC_agg_f, prot_silac_median)
prot_SILAC_agg_f$num_ratio <- 2^(log2(prot_SILAC_agg_f$num_ratio)-log2(prot_SILAC_agg_f$median))
prot_SILAC_agg_f <- prot_SILAC_agg_f[,-8]

check_median_step1 <- aggregate(data = prot_SILAC_agg_f, `num_ratio` ~ perc + sample + ratio, FUN = median, na.rm = TRUE) # check if median norm or scaling to ratio 1 worked

# continue downstream processing
prot_lf <- prot_SILAC_agg_f %>% dplyr::filter(`ratio` == "H/M") # select H/M SILAC ratios for quant comparison of IFNg treated cells / BSA control cells
prot_lf$log2FC <- log2(prot_lf$`num_ratio`)

# check distributions of log2 fold change values
ggplot(prot_lf, aes(x = as.factor(`sample`),y = `log2FC`, fill = as.factor(`sample`))) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.7) +
  geom_violin(show.legend = FALSE, alpha = 0.6) + xlab(NULL) + ylab("log2(SILAC H/M)") + 
  geom_boxplot(show.legend = FALSE, alpha = 0.6, width = 0.2) + theme_bw() + scale_x_discrete(limits = c("2h","4h","6h","9h","24h")) +
  theme(plot.title = element_text(color="black", size = 16, face= "bold"), axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(vjust = 0.5, color = "black", size = 14, angle = 0, hjust = 0.5),
        axis.text.y = element_text(color = "black", size = 14), axis.title.y = element_text(color = "black", size = 14, face = "bold"))


# calculate additional quant metrics
sd_l <- aggregate(data = prot_SILAC_agg_f, `num_ratio` ~ `perc` + `Genes` + `Protein.Ids` + `sample` + `ratio`, FUN = sd, na.rm = TRUE) # select corret ratio!
colnames(sd_l) <- c("perc","Gene_names","Protein.Group","sample","ratio","sd")

mean_l <- aggregate(data = prot_SILAC_agg_f, `num_ratio` ~ `perc` + `Genes` + `Protein.Ids` + `sample` + `ratio`, FUN = mean, na.rm = TRUE) # select corret ratio!
colnames(mean_l) <- c("perc","Gene_names","Protein.Group","sample","ratio","mean")

stats_l <- merge(mean_l, sd_l)
stats_l$cv <- (stats_l$sd / stats_l$mean)*100

# save for respective SILAC ratio
write_tsv(stats_l, "proteinGroups_stats_l_complete.tsv")
stats_l <- read_tsv(file = "proteinGroups_stats_l_complete.tsv")

# CV 
stats_lf <- stats_l #%>% dplyr::filter(`sample` == "AAAR1") # select indivdual sample for plot and plotgrid

p_meds <- aggregate(list(stats_lf$`cv`), by = list(stats_lf$`perc`, stats_lf$`sample`, stats_lf$`ratio`), median, na.rm = TRUE)
colnames(p_meds) <- c("perc","sample","ratio","cv")

#write_tsv(p_meds, "proteinGroups_median_CV.tsv")

middle <- ggplot(data = stats_lf %>% dplyr::filter(ratio == "H/M"), aes(y = `cv`, x = `sample`, fill = `sample`)) +
  geom_boxplot(show.legend = FALSE, alpha = 0.7, outlier.alpha = 0.1) + #scale_y_continuous(limits = c(0,1)) +
  geom_text(data = p_meds %>% dplyr::filter(ratio == "H/M"), aes(y = -10, x = `sample`, group = `sample`, label = round(`cv`, digits = 1)), col = "black", alpha = 0.9, size = 5, vjust = 0.5, 
            color = "black", position = position_dodge(width = 0.95), angle = 90) + scale_x_discrete(limits = c("2h","4h","6h","9h","24h")) +
  theme_bw() + #facet_grid(~perc) + #ggtitle("SILAC mix 1") +
  xlab("IFNg treatment") + ylab("coefficient of variation [%]") + ylim(c(-25,150)) + #scale_fill_manual(values = c("dodgerblue3","firebrick3")) +
  theme(plot.title = element_text(color="black", size = 16, face= "bold", vjust = 0.5), axis.title.x = element_text(color="black", size=14, face= "bold"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color="black", size = 14), axis.text.y = element_text(color="black", size=14), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(), axis.title.y = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 14, face = "bold"), legend.position = "right",
        legend.text = element_text(color="black", size=14), legend.title = element_text(color="black", size=14, face="bold"))

middle



prot_lf$`log2FC` <- log2(prot_lf$`num_ratio`)

proteins_wide <- prot_lf[,-c(4:8)] %>% spread(key = Experiment, value = `log2FC`)


## principle component analysis
df_lfc <- proteins_wide[,-c(2)]

numeric_only <- unlist(lapply(df_lfc, is.numeric)) # identify numerical columns in dataframe
df_LFC <- df_lfc %>% filter(rowSums(is.na(df_lfc[,numeric_only])) == 0)

# prepare df with numeric values (log2 fc valus) and protein ids as rownames
df_LFC <- df_LFC 
df_LFC <- df_LFC %>% remove_rownames %>% column_to_rownames(var = "Protein.Ids")

# meta data for PCA and tSNE
n_unique <- as.numeric(length(unique(prot_lf$sample)))
`treatment` = rep(c(unique(prot_lf$sample)), seq(1, 3, by = 1), length.out = 3*n_unique, each = 3) # multiply number of replicates with number of unique sample conditions

set.seed(123)
pca_data = prcomp(t(df_LFC[]))
pca_data_perc = round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
df_pca_data = data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2], sample = colnames(df_LFC), condition = `treatment`)
df_pca_data$complex <- stringr::str_extract(df_pca_data$condition, "^.{4}")


labs <- df_pca_data %>% group_by(condition) %>% filter(PC1 == max(PC1))

PCA <- ggplot(df_pca_data, aes(x = PC1, y = PC2, color = as.factor(`condition`))) + geom_point(size = 4, alpha = 0.7, show.legend = FALSE) + theme_bw() + 
  geom_encircle(size = 2, alpha = 0.6, expand = 0.05, show.legend = FALSE) + 
  geom_text_repel(data = labs, aes(x = PC1, y = PC2, label = `condition`), size = 4, show.legend = FALSE, color = "black") +
  #ggplot(df_pca_data, aes(x = PC1, y = PC2, color = as.factor(`condition`), shape = `group`)) + geom_point(size = 5, alpha = 0.7) + theme_bw() + geom_encircle(size = 2, alpha = 0.6, expand = 0.05) +
  labs(x=paste0("PC1 (",pca_data_perc[1],"%)"), y = paste0("PC2 (",pca_data_perc[2],"%)")) + labs(color = "target") +
  guides(color = guide_legend(title = "target gene", ncol = 4)) +
  #scale_color_manual(values = c("#DB72FB","#F8766D","#93AA00","#00BA38","#00B9E3")) +
  theme(plot.title = element_text(color="black", size=14, face= "bold", hjust=0.5), axis.title.x = element_text(color="black", size=14, hjust=0.5, face="bold"),
        axis.text.x = element_text(color="black", size=14), axis.text.y = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"), legend.text = element_text(color="black", size=14),
        legend.title = element_text(color="black", size=14, face="bold"), legend.position = "right",
        legend.background = element_blank()) #+ scale_color_manual(values = c("firebrick2", "mediumpurple3", "limegreen", "dodgerblue"))

PCA

plot_grid(quant_prot, middle, PCA, ncol = 3)

# save and import preursor count table
prec_count <- prec_count %>% spread(key = Experiment, value = `n`)

write_tsv(prec_count, "precursor_count.tsv")
prec_count <- read_tsv(file = "precursor_count.tsv")

# save and import wide protein group table
write_tsv(proteins_wide, "proteinGroups_wide_df.tsv")
proteins_wide <- read_tsv(file = "proteinGroups_wide_df.tsv")
proteins_wide <- dplyr::filter(proteins_wide, !grepl(";", Protein.Ids)) # remove proteins groups with non unique identifier for differential expression analysis

# count missing values
numeric_only <- unlist(lapply(proteins_wide, is.numeric)) # identify numerical columns in dataframe
pg_ratio_mat = as.matrix(proteins_wide[,numeric_only])  ### select samples of interest !!!!
rownames(pg_ratio_mat) <- proteins_wide$Protein.Ids
na_pg_ratios <- sum(is.na(pg_ratio_mat))
valid_pg_ratios <- sum(!is.na(pg_ratio_mat))
total <- na_pg_ratios + valid_pg_ratios
perc_NA <- (na_pg_ratios/total)*100



# differential expression analysis
sample_select <- "2h" # select time point to subset for differential expression testing

# diff abundance test
df_subset <- proteins_wide %>% dplyr::select(`Protein.Ids`, starts_with(paste0(sample_select,"_R")))

df_subset[df_subset == 0] <- NA

numeric_only <- unlist(lapply(df_subset, is.numeric)) # identify numerical columns in dataframe
df_subset <- df_subset %>% filter(rowSums(is.na(df_subset[,numeric_only])) <= 1)  %>% distinct(`Protein.Ids`, .keep_all = TRUE) # 2 in 3 replicates

protein_df <- df_subset
numeric_only <- unlist(lapply(protein_df, is.numeric)) # identify numerical columns in dataframe

protein_matrix = as.matrix(protein_df[,numeric_only])  ### select samples of interest !!!!
rownames(protein_matrix) = protein_df$Protein.Ids

protein_matrix <- scale(protein_matrix, scale = FALSE, center = TRUE) #  median normalize log2 ratio values to 0

# without contrast fir
fit2_ref = lmFit(protein_matrix)
fit3_ref <- eBayes(fit2_ref, robust = TRUE)


# multiscatter plot for log2 SILAC ratios - low correlation with distribution of values arund 0 indicates small effect of treatment
multiscatter_function <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = proteins_wide, mapping = mapping) + geom_point(alpha = 0.1, colour = "black", size = 2) +
    geom_abline(mapping = NULL, data = NULL, ..., slope = 1, intercept = 0, na.rm = FALSE, show.legend = NA) +
    geom_smooth(method = method, color = "red", ...)
 p }

ggpairs(data = protein_df, columns = c(2:4), lower = list(continuous = wrap(multiscatter_function, method = "lm")),
       diag = list(continuous = wrap("barDiag", colour = "black", bins = 75)),
       upper = list(continuous = wrap("cor", method="pearson", size = 5))) + theme_bw() + theme(strip.text = element_text(size = 12, face = "bold"))


# create table of precurosr counts to include for DEqMS
peptide_count_subset <- prec_count %>% dplyr::select(`Protein.Ids`, contains(paste0(sample_select,"_R"))) %>% distinct(`Protein.Ids`, .keep_all = TRUE)
sum(duplicated(peptide_count_subset$Protein.Ids))
tool_df <- as.data.frame(protein_df[,c("Protein.Ids")])
tool_df$junk <- rep("X")
colnames(tool_df) <- c("Protein.Ids","junk")

peptide_count_subset <- full_join(tool_df, peptide_count_subset)
peptide_count_subset <- peptide_count_subset[,-2] %>% distinct(`Protein.Ids`, .keep_all = TRUE)
peptide_count <- peptide_count_subset

numeric_only <- unlist(lapply(peptide_count, is.numeric)) # identify numerical columns in dataframe

peptide_count <- as.data.frame(peptide_count)
rownames(peptide_count) <- peptide_count$Protein.Ids
peptide_count <- peptide_count[,numeric_only]

pep.count.table = data.frame(count = rowMins(as.matrix(peptide_count)), row.names = rownames(peptide_count))
pep.count.table[is.na(pep.count.table)] <- 0
pep.count.table$count = pep.count.table$count + 1

fit3_ref$count = pep.count.table[rownames(fit3_ref$coefficients),"count"]

#check the values in the vector fit3$count
#if min(fit3$count) return NA or 0, you should troubleshoot the error first
min(fit3_ref$count)
max(fit3_ref$count)

fit4_ref = spectraCounteBayes(fit3_ref, fit.method = "spline")

# extract DEqMS outputs
DEqMS.results = outputResult(fit4_ref, coef_col = 1)
head(DEqMS.results)

DEqMS_result <- DEqMS.results
DEqMS_result <- tibble::rownames_to_column(DEqMS_result, "Protein.Ids")

DEqMS_result <- merge(DEqMS_result, proteins_wide[,c("Genes","Protein.Ids")], by = c("Protein.Ids"))
sum(duplicated(DEqMS_result$Protein.Ids))
DEqMS_result <- DEqMS_result %>% distinct(`Protein.Ids`,`logFC`,`sca.adj.pval`, .keep_all = TRUE)

# add average expression in control to differential expression output table
prot_lfq <- read_tsv("iq_pg_proc_MS2.tsv")
df_ctrl_int <- prot_lfq %>% dplyr::select(`Protein.Ids`, starts_with(paste0(sample_select,"_M"))) #  control condition in this case is in intermediate SILAC channel
numeric_only <- unlist(lapply(df_ctrl_int, is.numeric)) # identify numerical columns in dataframe
df_ctrl_int$AveExpr <- rowMeans(df_ctrl_int[,numeric_only], na.rm = TRUE)

DEqMS_result <- merge(DEqMS_result[,-3], df_ctrl_int[,c("Protein.Ids","AveExpr")], by = c("Protein.Ids"))

# save table(s) and repeat step for selected timepoint

#write_tsv(DEqMS_result, "DEqMS_IFNg_NP_2h.tsv")
DEqMS_result <- read_tsv(file = "DEqMS_IFNg_NP_2h.tsv")
#write_tsv(DEqMS_result, "DEqMS_IFNg_NP_4h.tsv")
DEqMS_result <- read_tsv(file = "DEqMS_IFNg_NP_4h.tsv")
#write_tsv(DEqMS_result, "DEqMS_IFNg_NP_6h.tsv")
DEqMS_result <- read_tsv(file = "DEqMS_IFNg_NP_6h.tsv")
#write_tsv(DEqMS_result, "DEqMS_IFNg_NP_9h.tsv")
DEqMS_result <- read_tsv(file = "DEqMS_IFNg_NP_9h.tsv")
#write_tsv(DEqMS_result, "DEqMS_IFNg_NP_24h.tsv")
DEqMS_result <- read_tsv(file = "DEqMS_IFNg_NP_24h.tsv")

# 
sample_select <- "2h" # select time point to subset for differential expression testing

# plot
volcano_DEqMS <- ggplot(data = DEqMS_result, aes(x = logFC, y =  -log10(`sca.adj.pval`))) + 
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.7, size = 1, color = "darkgrey") +
  geom_point(alpha = 0.1, size = 2) + #facet_wrap(~treatment, ncol = 4) +
  geom_point(data = DEqMS_result %>% filter(`logFC` > 0.585 & `sca.adj.pval` < 0.05), color = "firebrick3", size = 2, alpha = 0.3) +
  geom_point(data = DEqMS_result %>% filter(`logFC` < -0.585 & `sca.adj.pval` < 0.05), color = "dodgerblue3", size = 2, alpha = 0.3) +
  geom_text_repel(color = "black", data = DEqMS_result %>% filter(`sca.adj.pval` < 0.01 & abs(`logFC`) > 1.0), mapping = aes(label = `Genes`), size = 4) + 
  #geom_point(data = DEqMS_result %>% filter(`gene` == "P60842;Q14240"), color = "green", size = 3, alpha = 0.7) +
  #geom_label_repel(color = "black", data = DEqMS_result %>% filter(`gene` == "P60842;Q14240"), mapping = aes(label = `Genes`), size = 4) + 
  ggtitle(paste0(sample_select," - IFNg")) + 
  xlab(paste0("log2(IFNg / control)")) + ylab("-log10(adj. p-value)") + theme_bw() + 
  xlim(c(-3.4,5.4)) + ylim(c(0,4.25)) +
  theme(plot.title = element_text(color="black", size=14, face="bold"), 
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_blank(), axis.text.x = element_text(color="black", size=14), axis.text.y = element_blank(), strip.text = element_blank(),
        #axis.text.x = element_text(color="black", size=14), axis.text.y = element_text(color="black", size=14), strip.text = element_blank(), axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title = element_text(color="black", size=14, face="bold"), 
        legend.position = "bottom", legend.text = element_text(size = 14))

x_dens <- axis_canvas(volcano_DEqMS, axis = "x") + geom_density(data = DEqMS_result, aes(x = `logFC`), alpha = 0.7, fill = "black") + geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.7, size = 1, color = "darkgrey") + scale_fill_manual(values = c("#009E73","#000000"))#+ geom_vline(xintercept = median(DEqMS_result$logFC), alpha = 0.5, size = 0.5, color = "black")

p1 <- insert_xaxis_grob(volcano_DEqMS, x_dens, grid::unit(0.15, "null"), position = "top")
ggdraw(p1)

# name volcano plot of the 5 different timepoints and plot in a panel
h1 <- ggdraw(p1)
h2 <- ggdraw(p1)
h3 <- ggdraw(p1)
h4 <- ggdraw(p1)
h5 <- ggdraw(p1)

volcano_panel <- plot_grid(h1,h2,h3,h4,h5, ncol = 5, nrow = 1, align = "hv", labels = c('2h','4h','6h','9h','24h'), label_size = 14, label_y =  1.1,
                           rel_widths = c(0.22,0.19,0.19,0.19,0.19)) # set x axis values in density and volcano
volcano_panel


# full panel
plot_grid(panel_1, volcano_panel, eGO_dot, ncol = 1, align = "v")




### enrichment analysis
# import tables with different identifyerse.g Uniprot, Ensembl, Entrez and gene symbols for human proteome
uniprot_ID <- read_tsv(file = "Hs_uniprot_IDs.tsv")
uniprot_ID$ENTREZID <- as.character(uniprot_ID$ENTREZID)

anno_tbl <- read_tsv(file = "Hs_annotation_df.tsv")
anno_tbl$ENTREZID <- as.character(anno_tbl$ENTREZID)

# RNAseq data - 24h IFNg
top_500_RNA_up <- read_tsv(file = "Hela_24h_IFNg_top_500_RNA_up.tsv") # top 500 unpregulated RNAs from "GSE150196_RNA-seq_DESeq2_priming_vs_naive.tab"

# ChIPseq STAT1 30 min IFNg
top_500_ChIP <- read_tsv(file = "Hela_30min_IFNg_top_500_ChIP.tsv") #  top 500 ChIPseq targts of STAT1 in Hela cells with 30 min IFNg stimulation - ENCODE data ENCFF039MZH


## retrieve MolSign Database genelists ###
# create term2gene dataframe , category = "H" == hallmark gene sets , "C2" == curated gene sets , "C3" == motif genesets
HALLMARK <- msigdbr(species = "Homo sapiens", category = c("H")) %>% dplyr::select(gs_name, entrez_gene) # hallmark genesets
C2 <- msigdbr(species = "Homo sapiens", category = c("C2")) %>% dplyr::select(gs_name, entrez_gene) # curated gene sets
C5 <- msigdbr(species = "Homo sapiens", category = c("C5")) %>% dplyr::select(gs_name, entrez_gene) # ontology gene sets

m_t2g <- rbind(HALLMARK, top_500_RNA_up[,c("gs_name", "entrez_gene")], top_500_ChIP[,c("gs_name", "entrez_gene")])
#m_t2g <- rbind(HALLMARK, C2, C5)
head(m_t2g)

up_IDs <- uniprot_ID
colnames(up_IDs)[colnames(up_IDs) == "ENTREZID"] <- "entrez_gene"
m_t2g_2 <- merge(m_t2g, up_IDs, all.x = TRUE)


# feed in your differential expression output table
backgroundgenes.df <- DEqMS_result # load the DEqMS output table for the selectied time point !!!
backgroundgenes.df <- backgroundgenes.df #%>% mutate(Protein.Ids = str_split(Protein.Ids,";", simplify = TRUE)[,1]) # simplify protein Ids to first uniprot Id for GSEA - not needed here since we only use unique protein groups
colnames(backgroundgenes.df)[colnames(backgroundgenes.df)=="Protein.Ids"] <- "UNIPROT" # rename to merge with referenc tables and add Entrez Ids
#backgroundgenes.df <- backgroundgenes.df %>% distinct(`Protein.Ids`,`logFC`,`sca.adj.pval`, .keep_all = TRUE)

genelist_v3 <- merge(backgroundgenes.df, uniprot_ID, by="UNIPROT", all.x = TRUE)
GSEA_genelist_v3 <- genelist_v3$logFC
names(GSEA_genelist_v3) = as.character(genelist_v3$UNIPROT)
GSEA_genelist_v3 = sort(GSEA_genelist_v3, decreasing = TRUE)

register(DoparParam()) # avoid errors of parellelization via BiocParallel


### ORA (overrepresentation enrichment analysis)
DOWNreg <- genelist_v3 %>% dplyr::filter(`sca.adj.pval` < 0.05 & `logFC` < -0.585)
UPreg <- genelist_v3 %>% dplyr::filter(`sca.adj.pval` < 0.05 & `logFC` > 0.585)

## MSigDB ORA
# Downregulated proteins
# only single protein for 2h timepoint gives error
em_down <- enricher(DOWNreg$UNIPROT, TERM2GENE = m_t2g_2[,2:3], pAdjustMethod = "BH", universe = genelist_v3$UNIPROT, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.05)
enricher_result_down <- em_down@result
enricher_result_down <- mutate(enricher_result_down, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
enricher_result_down <- mutate(enricher_result_down, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
enricher_result_down$direction <- rep("down")
  
# Upregulated proteins
em_up <- enricher(UPreg$UNIPROT, TERM2GENE = m_t2g_2[,2:3], pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = genelist_v3$UNIPROT, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.05)
enricher_result_up <- em_up@result
enricher_result_up <- mutate(enricher_result_up, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
enricher_result_up <- mutate(enricher_result_up, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
enricher_result_up$direction <- rep("up")

enricher_result <- rbind(enricher_result_up, enricher_result_down) # 
#enricher_result <- rbind(enricher_result_up) # no multiple downregulated protein for time point 2h 

# save results of ORA after going through the different timepoints -  repeat previous steps!
write_tsv(enricher_result, "HALLMARK_comp_ORA_2h.tsv")
write_tsv(enricher_result, "HALLMARK_comp_ORA_4h.tsv")
write_tsv(enricher_result, "HALLMARK_comp_ORA_6h.tsv")
write_tsv(enricher_result, "HALLMARK_comp_ORA_9h.tsv")
write_tsv(enricher_result, "HALLMARK_comp_ORA_24h.tsv")

## GO over representation enrichment analysis (slower)
ego_up <- enrichGO(gene = UPreg$UNIPROT, keyType = "UNIPROT", universe = genelist_v3$UNIPROT, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
ego_result_up <- ego_up@result
ego_result_up$direction <- rep("up")

ego_down <- enrichGO(gene = DOWNreg$UNIPROT, keyType = "UNIPROT", universe = genelist_v3$UNIPROT, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
ego_result_down <- ego_down@result
ego_result_down$direction <- rep("down")

comb_eGO <- rbind(ego_result_up, ego_result_down)
comb_eGO$name <- paste0(comb_eGO$Description, "-", comb_eGO$ID)

# save results of GO ORA after going through the different timepoints -  repeat previous steps!
write_tsv(comb_eGO, "eGO_ORA_2h.tsv")
write_tsv(comb_eGO, "eGO_ORA_4h.tsv")
write_tsv(comb_eGO, "eGO_ORA_6h.tsv")
write_tsv(comb_eGO, "eGO_ORA_9h.tsv")
write_tsv(comb_eGO, "eGO_ORA_24h.tsv")

## GSEA
em2 <- GSEA(GSEA_genelist_v3, TERM2GENE = m_t2g_2[,c("gs_name","UNIPROT")], pvalueCutoff = 0.20, minGSSize = 10, maxGSSize = 500, eps = 0, exponent = 1, seed = T, pAdjustMethod = "BH")
GSEA_result_molSig <- em2@result

gseaplot(em2, geneSetID = "HALLMARK_INTERFERON_GAMMA_RESPONSE", by = "runningScore", title = "HALLMARK_INTERFERON_GAMMA_RESPONSE")

# GO-GSEA
ego3 <- gseGO(geneList = GSEA_genelist_v3, OrgDb = org.Hs.eg.db, keyType = "UNIPROT",  ont = "ALL", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, eps = 0, seed = T)
GO_GSEA_res <- ego3@result

gseaplot(ego3, geneSetID = "GO:0034341", by = "runningScore", title = "response to interferon-gamma - GO:0034341")


## highlight temporal changes upon IFNg stimulation
DEqMS_result_2h <- read_tsv(file = "DEqMS_IFNg_NP_2h.tsv")
DEqMS_result_2h$time <- rep(2)
DEqMS_result_4h <- read_tsv(file = "DEqMS_IFNg_NP_4h.tsv")
DEqMS_result_4h$time <- rep(4)
DEqMS_result_6h <- read_tsv(file = "DEqMS_IFNg_NP_6h.tsv")
DEqMS_result_6h$time <- rep(6)
DEqMS_result_9h <- read_tsv(file = "DEqMS_IFNg_NP_9h.tsv")
DEqMS_result_9h$time <- rep(9)
DEqMS_result_24h <- read_tsv(file = "DEqMS_IFNg_NP_24h.tsv")
DEqMS_result_24h$time <- rep(24)

prot_SILAC_agg <- read_tsv(file = "aggregated_SILAC_ratios.tsv")
prot_SILAC_agg <- prot_SILAC_agg[!is.na(prot_SILAC_agg$num_ratio),] # remove 0 values for norm SILAC ratios
prot_SILAC_agg$n <- rep(1)
prot_SILAC_agg_f <- prot_SILAC_agg %>% dplyr::filter(`ratio` == "H/M") #%>% arrange(sample) %>% mutate(sample = factor(sample, levels = c("2h", "4h", "6h", "9h", "24h")))
prot_SILAC_agg_f2 <- prot_SILAC_agg_f
prot_SILAC_agg_f2$time = as.numeric(substr(prot_SILAC_agg_f2$sample,1,nchar(prot_SILAC_agg_f2$sample)-1))

long_DEqMS <- rbind(DEqMS_result_2h, DEqMS_result_4h, DEqMS_result_6h, DEqMS_result_9h, DEqMS_result_24h)
long_DEqMS$group <- rep("all proteins")

# reference list of IFNg targets
IFNg_ref_set <- m_t2g_2[,2:3] %>% dplyr::filter(`gs_name` == "GOBP_RESPONSE_TO_INTERFERON_GAMMA"|`gs_name` == "HALLMARK_INTERFERON_GAMMA_RESPONSE")
colnames(IFNg_ref_set) <- c("gs_name","Protein.Ids")

# alternative selection of protein traces
grouped_prot <- long_DEqMS %>% dplyr::filter(`sca.adj.pval` < 0.05 & abs(`logFC`) > 1)

grouped_prot_l <- inner_join(long_DEqMS, unique(grouped_prot[,c("Protein.Ids","Genes")]))
POI_deq <- grouped_prot_l %>% dplyr::filter(!grepl(";", Protein.Ids)) # remove proteins groups with non unique identifier
POI_subset <- merge(prot_SILAC_agg_f2, grouped_prot[,c("Protein.Ids","Genes")])
POI_subset <- dplyr::filter(POI_subset, !grepl(";", Protein.Ids)) # remove proteins groups with non unique identifier

sd_set <- aggregate(data = POI_subset, `num_ratio` ~ `perc` + `Genes` + `Protein.Ids` + `sample` + `ratio`, FUN = sd, na.rm = TRUE) # select corret ratio!
colnames(sd_set) <- c("perc","Genes","Protein.Ids","sample","ratio","sd")
mean_set <- aggregate(data = POI_subset, `num_ratio` ~ `perc` + `Genes` + `Protein.Ids` + `sample` + `ratio`, FUN = mean, na.rm = TRUE) # select corret ratio!
colnames(mean_set) <- c("perc","Genes","Protein.Ids","sample","ratio","mean")
stats_set <- merge(mean_set, sd_set)
stats_set$cv <- (stats_set$sd / stats_set$mean)*100
stats_set <- aggregate(data = stats_set, `cv` ~ `perc` + `Genes` + `Protein.Ids`, FUN = max, na.rm = TRUE) # get max CV values per protein groups at all timepoints
stats_set <- stats_set %>% dplyr::filter(cv < 20) # remove protein groups with large cV values at any time point

POI_subset <- merge(POI_subset, stats_set[,c("Protein.Ids","Genes")])
POI_subset <- unique(POI_subset)
POI_deq <- merge(POI_deq, stats_set[,c("Protein.Ids","Genes")])
POI_deq <- unique(POI_deq)

POI_subset_IFNg <- inner_join(POI_subset, IFNg_ref_set)
POI_subset_IFNg <- POI_subset_IFNg[,1:9]
POI_subset_IFNg <- unique(POI_subset_IFNg)
POI_subset_IFNg$group <- rep("annotated IFNg targets*")
POI_deq_IFNg <- inner_join(POI_deq, IFNg_ref_set)
POI_deq_IFNg <- POI_deq_IFNg[,1:15]
POI_deq_IFNg <- unique(POI_deq_IFNg)
POI_deq_IFNg$group <- rep("annotated IFNg targets*")

POI_subset_rest <- anti_join(POI_subset, IFNg_ref_set)
POI_subset_rest <- POI_subset_rest[,1:9]
POI_subset_rest <- unique(POI_subset_rest)
POI_subset_rest$group <- rep("non-annotated IFNg targets*")
POI_deq_rest <- anti_join(POI_deq, IFNg_ref_set)
POI_deq_rest <- POI_deq_rest[,1:15]
POI_deq_rest <- unique(POI_deq_rest)
POI_deq_rest$group <- rep("non-annotated IFNg targets*")

POI_subset <- rbind(POI_subset_IFNg, POI_subset_rest)
POI_deq <- rbind(POI_deq_IFNg, POI_deq_rest)

POI_deq <- POI_deq %>% add_column(`is_significant` = if_else(POI_deq$sca.adj.pval < 0.05 & abs(POI_deq$logFC) > 0.585, TRUE, FALSE))
POI_deq$time_h <- paste0(POI_deq$time,"h")

# include subgroups for heatmap depending on significance at tme points
early <- POI_deq %>% dplyr::filter(abs(logFC) > 0.582 & sca.adj.pval < 0.05 & time <= 2)
early$subgroup <- rep("early")
early <- inner_join(POI_deq, early[,c("Protein.Ids","subgroup")])

intermediate <- POI_deq %>% dplyr::filter(abs(logFC) > 0.582 & sca.adj.pval < 0.05 & time >= 4 & time <= 9)
intermediate$subgroup <- rep("intermediate")
intermediate <- inner_join(POI_deq, intermediate[,c("Protein.Ids","subgroup")])
intermediate <- anti_join(intermediate, early[,c("Protein.Ids","Genes")])

late <- POI_deq %>% dplyr::filter(abs(logFC) > 0.582 & sca.adj.pval < 0.05 & time >= 24)
late$subgroup <- rep("late")
late <- inner_join(POI_deq, late[,c("Protein.Ids","subgroup")])
late <- anti_join(late, early[,c("Protein.Ids","Genes")]) %>% anti_join(intermediate[,c("Protein.Ids","Genes")])

POI_deq <- rbind(early, intermediate, late)

# create reference gene list for search on interferome.org
write_tsv(unique(POI_deq[,c("Protein.Ids","Genes","group","subgroup")]), "reduced_POI_deq_list.tsv")

# import interferome data
interferome <- read_tsv(file = "interferome_DataSearchResults_complete.txt")
interferome2 <- read_tsv(file = "interferome_DataSearchResults_complete_synonyms_of_missing.txt")
interferome <- rbind(interferome, interferome2)
interferome_a <- interferome[,c(1:5)]
interferome_a$n_datasets_interferome <- rep(1)
n_ifn <- aggregate(data = interferome_a, n_datasets_interferome ~ Genes, FUN = sum)
n_ifn <- unique(n_ifn)

n_ifn[n_ifn$Genes == "WARS", "Genes"] <- "WARS1" # need to manually change Gene names of several proteins..... thanks interferome.org
n_ifn[n_ifn$Genes == "ADCK4", "Genes"] <- "COQ8B" # need to manually change Gene names of several proteins..... thanks interferome.org
n_ifn[n_ifn$Genes == "FAM105A", "Genes"] <- "OTULINL" # need to manually change Gene names of several proteins..... thanks interferome.org
n_ifn[n_ifn$Genes == "GRAMD3", "Genes"] <- "GRAMD2B" # need to manually change Gene names of several proteins..... thanks interferome.org
n_ifn[n_ifn$Genes == "PPAP2B", "Genes"] <- "PLPP3" # need to manually change Gene names of several proteins..... thanks interferome.org
n_ifn[n_ifn$Genes == "PVRL2", "Genes"] <- "NECTIN2" # need to manually change Gene names of several proteins..... thanks interferome.org
n_ifn[n_ifn$Genes == "TMEM173", "Genes"] <- "STING1" # need to manually change Gene names of several proteins..... thanks interferome.org


POI_deq <- full_join(POI_deq, n_ifn)
POI_deq <- unique(POI_deq)

# add number of datapoints per protein group
POI_deq$n_datapoints <- rep(1)
n_datapoints <- aggregate(data = POI_deq, n_datapoints ~ Protein.Ids + Genes, FUN = sum)

#POI_deq <- left_join(POI_deq[,-c(19,21)], n_datapoints)
POI_deq <- left_join(POI_deq[,-c(20)], n_datapoints)
POI_deq$group2 <- rep("differentially expressed proteins")
POI_deq <- POI_deq %>% add_column(`GO or Hallmark` = if_else(POI_deq$group == "annotated IFNg targets*", TRUE, FALSE))

# save and import table
write_tsv(POI_deq, "selected_diffregPROT_IFNg_heatmap.tsv")
POI_deq <- read_tsv(file = "selected_diffregPROT_IFNg_heatmap.tsv")

# plot heatmap
heatmap <- ggplot(data = POI_deq, aes(x = as.factor(`time_h`), y = reorder(`Genes`, `logFC`))) + 
  geom_tile(mapping = aes( fill = `logFC`)) + 
  geom_tile(data = POI_deq %>% dplyr::filter(`is_significant` == TRUE), mapping = aes( fill = `logFC`), color = "black", linewidth = 0.1) + 
  #geom_text(POI_deq, mapping = aes(x = as.factor(`time_h`), y = reorder(`Genes`, `logFC`), label = round(`logFC`, digits = 2)), color = "black", size = 3) +  
  scale_y_discrete(position = "right") + scale_x_discrete(limits = c("2h","4h","6h","9h","24h")) + labs(fill = 'log2(FC)') + xlab("IFNg treatment") +
  scale_fill_gradient2(low = "dodgerblue3", mid = "white", high = "firebrick3", limits = c(-2.5,2.5), na.value = "firebrick3") + 
  facet_nested(vars(`group2`,`subgroup`), strip = strip_nested(bleed = TRUE, size = "variable", by_layer_y = FALSE, 
                                                              background_y = elem_list_rect(fill = c("lightgray","gold","limegreen","orchid"))),
               scales = "free", switch = "y") + force_panelsizes(rows = c(6,42,55)) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 10), axis.title = element_blank(), panel.grid = element_blank(),
        strip.text = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14), legend.position = "bottom", legend.spacing.y = unit(0.5, "mm"),
        legend.title = element_text(size = 14, face = "bold"), axis.title.x = element_text(size = 14, face = "bold"), axis.title.y = element_blank())

#heatmap

n_ifn_ome <- ggplot(data = POI_deq, aes(fill = `subgroup`, x = (`n_datasets_interferome`/`n_datapoints`), y = reorder(`Genes`, `logFC`))) + 
  geom_col(show.legend = TRUE) + labs(fill = 'sub group') + xlab("reports on interferome.org") + scale_x_continuous(limits = c(-5,155)) +
  geom_point(data = POI_deq %>% dplyr::filter(`GO or Hallmark` == TRUE), aes(x = `n_datasets_interferome`+15, y = reorder(`Genes`, `logFC`)), size = 1.5, shape = 8, alpha = 0.8, show.legend = TRUE) + 
  geom_text(POI_deq, mapping = aes(x = `n_datasets_interferome`-5, y = reorder(`Genes`, `logFC`), label = `n_datasets_interferome`), color = "black", size = 3.5) +  
  facet_nested(vars(`group2`,`subgroup`), strip = strip_nested(bleed = TRUE, size = "variable", by_layer_y = FALSE, 
                                                               background_y = elem_list_rect(fill = c("lightgray","gold","limegreen","orchid"))),
               scales = "free") + force_panelsizes(rows = c(6,42,55)) + scale_fill_manual(values = c("gold","limegreen","orchid")) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_blank(), axis.title.y = element_blank(), panel.grid = element_blank(), legend.spacing.y = unit(0.5, "mm"),
        strip.text = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14), legend.position = "bottom", legend.box="vertical", legend.margin=margin(),
        legend.title = element_text(size = 14, face = "bold"), axis.title.x = element_text(size = 14, face = "bold"))

n_ifn_ome <- n_ifn_ome+guides(fill = guide_legend(nrow = 3, byrow=TRUE))
#n_ifn_ome

plot_grid(heatmap, n_ifn_ome, ncol = 2, rel_widths = c(0.7,0.3), align = "h")

# create dotplot for overrepresentation enrichment results
# import GO results
GSEA_1 <- read_tsv(file = "eGO_ORA_2h.tsv") # does not exsist since there no significant enrichments at2h time point (with GO-terms)
GSEA_1$time <- rep("2h", nrow(GSEA_1))
GSEA_1$method <- rep("nasc. prot", nrow(GSEA_1))

GSEA_2 <- read_tsv(file = "eGO_ORA_4h.tsv")
GSEA_2$time <- rep("4h", nrow(GSEA_2))
GSEA_2$method <- rep("nasc. prot", nrow(GSEA_2))

GSEA_3 <- read_tsv(file = "eGO_ORA_6h.tsv")
GSEA_3$time <- rep("6h", nrow(GSEA_3))
GSEA_3$method <- rep("nasc. prot", nrow(GSEA_3))

GSEA_4 <- read_tsv(file = "eGO_ORA_9h.tsv")
GSEA_4$time <- rep("9h", nrow(GSEA_4))
GSEA_4$method <- rep("nasc. prot", nrow(GSEA_4))

GSEA_5 <- read_tsv(file = "eGO_ORA_24h.tsv")
GSEA_5$time <- rep("24h", nrow(GSEA_5))
GSEA_5$method <- rep("nasc. prot", nrow(GSEA_5))


# filter and combine df
comb_GSEA <- rbind(GSEA_2, GSEA_3, GSEA_4, GSEA_5)
comb_GSEA <- comb_GSEA %>% filter(`qvalue` < 0.05)
comb_GSEA <- mutate(comb_GSEA, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
comb_GSEA <- mutate(comb_GSEA, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
#comb_GSEA <- comb_GSEA %>% filter(`qvalues` < 0.05)
comb_GSEA$name <- paste0(comb_GSEA$Description, " - ", comb_GSEA$ID)

# save GO ORA table
write_tsv(comb_GSEA, "combined_eGO_results.tsv")
comb_GSEA <- read_tsv(file = "combined_eGO_results.tsv")

# repeat with Hallmark ORA results
GSEA_1 <- read_tsv(file = "HALLMARK_comp_ORA_2h.tsv")
GSEA_1$time <- rep("2h", nrow(GSEA_1))
GSEA_1$method <- rep("nasc. prot", nrow(GSEA_1))

GSEA_2 <- read_tsv(file = "HALLMARK_comp_ORA_4h.tsv")
GSEA_2$time <- rep("4h", nrow(GSEA_2))
GSEA_2$method <- rep("nasc. prot", nrow(GSEA_2))

GSEA_3 <- read_tsv(file = "HALLMARK_comp_ORA_6h.tsv")
GSEA_3$time <- rep("6h", nrow(GSEA_3))
GSEA_3$method <- rep("nasc. prot", nrow(GSEA_3))

GSEA_4 <- read_tsv(file = "HALLMARK_comp_ORA_9h.tsv")
GSEA_4$time <- rep("9h", nrow(GSEA_4))
GSEA_4$method <- rep("nasc. prot", nrow(GSEA_4))

GSEA_5 <- read_tsv(file = "HALLMARK_comp_ORA_24h.tsv")
GSEA_5$time <- rep("24h", nrow(GSEA_5))
GSEA_5$method <- rep("nasc. prot", nrow(GSEA_5))


# filter and combine df
comb_GSEA2 <- rbind(GSEA_1, GSEA_2, GSEA_3, GSEA_4, GSEA_5)
comb_GSEA2 <- comb_GSEA2 %>% filter(`qvalue` < 0.05)
comb_GSEA2$name <- comb_GSEA2$Description
comb_GSEA2$ONTOLOGY <- rep(NA)
comb_GSEA2$name[comb_GSEA2$name == "HALLMARK_APOPTOSIS"] <- "Hallmark apoptosis - M5902"
comb_GSEA2$name[comb_GSEA2$name == "HALLMARK_INTERFERON_GAMMA_RESPONSE"] <- "Hallmark Interferon gamma response - M5913"
comb_GSEA2$name[comb_GSEA2$name == "HALLMARK_ALLOGRAFT_REJECTION"] <- "Hallmark allograft rejection - M5950"

# combine separate gene set enrichment datasets
comb_GSEA <- rbind(comb_GSEA, comb_GSEA2)

# save full ORA result table
write_tsv(comb_GSEA, "full_ORA_enrichment.tsv")
comb_GSEA <- read_tsv(file = "full_ORA_enrichment.tsv")

# manually selected genesets for dotplot
GSEA_f <- rbind(#comb_GSEA %>% filter(`Description` == "GO_UNFOLDED_PROTEIN_BINDING"),
  comb_GSEA %>% filter(`Description` == "innate immune response"),
  #comb_GSEA %>% filter(`Description` == "defense response to virus"),
  #comb_GSEA %>% filter(`Description` == "type I interferon signaling pathway"),
  comb_GSEA %>% filter(`Description` == "interferon-gamma-mediated signaling pathway"),
  comb_GSEA %>% filter(`Description` == "response to interferon-gamma"),
  #comb_GSEA %>% filter(`Description` == "cellular response to interferon-gamma"),
  #comb_GSEA %>% filter(`Description` == "antigen processing and presentation of peptide antigen via MHC class I"),
  comb_GSEA %>% filter(`Description` == "antigen processing and presentation of endogenous antigen"),
  #comb_GSEA %>% filter(`Description` == "extracellular matrix structural constituent"),
  #comb_GSEA %>% filter(`Description` == "rRNA processing"),
  comb_GSEA %>% filter(`Description` == "top 500 upreg. RNA - GSE150196"),
  comb_GSEA %>% filter(`Description` == "top 500 STAT1 ChIPseq peaks - ENCFF039MZH"),
  comb_GSEA %>% filter(`Description` == "HALLMARK_APOPTOSIS"),
  #comb_GSEA %>% filter(`Description` == "HALLMARK_ALLOGRAFT_REJECTION"),
  comb_GSEA %>% filter(`Description` == "HALLMARK_INTERFERON_GAMMA_RESPONSE"))
#comb_GSEA %>% filter(`Description` == "interferon-gamma-mediated signaling pathway"))


# dot plot of ORA
# with grid
legend_title1 <- c("fold enrichment")
legend_title2 <- c("-log10(q-value)")
#legend_title3 <- c("ontology")

eGO_dot <- ggplot() + 
  #geom_point(GSEA_f, mapping = aes(x = `name`, y = `time`, size = `FoldEnrichment`, shape = `ONTOLOGY`), stat = 'identity', color = "black", stroke = 2) +
  geom_point(GSEA_f, mapping = aes(x = `name`, y = `time`, size = `FoldEnrichment`, color = -log10(`qvalue`))) +
  geom_point(GSEA_f, mapping = aes(x = `name`, y = `time`, size = `FoldEnrichment`), shape = 21, color = "black", stroke = 2) +
  ylab("IFN-g treatment") + xlab("") + #scale_shape_discrete(legend_title3) +
  scale_radius(legend_title1, range = c(3,9), breaks = c(10,30)) + theme_bw() + #facet_grid(~method, scales = "free") +
  scale_color_distiller(legend_title2, palette = "OrRd", direction = 1) +
  scale_y_discrete(limits = c("2h","4h","6h","9h","24h")) + coord_flip() + #facet_grid(~group) +
  theme(plot.title = element_text(color="black", size=14, face= "bold", hjust=0.5), axis.title.x = element_text(color="black", size=14, hjust=0.5, face="bold"), axis.text.x = element_text(color="black", size=14),
        axis.text.y = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14), legend.text = element_text(color="black", size=14),
        legend.title = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 14, face = "bold"))

eGO_dot


# correlation scatter
DEqMS_result <- read_tsv(file = "DEqMS_IFNg_NP_24h.tsv")
DEqMS_result <- DEqMS_result %>% mutate(`Protein.Ids` = strsplit(as.character(Protein.Ids), ";")) %>% unnest(Protein.Ids)

# RNA-seq
DEseq2_24h <- read_tsv(file = "DEseq2_24h_IFNg.tsv")
DEseq2_24h <- rename(DEseq2_24h, c('UNIPROT' = 'Protein.Ids'))

STAT1_ChIPseq_targets <- read_tsv(file = "STAT1_ChIPseq_targets.tsv")
no_dup_ChIP <- aggregate(data = STAT1_ChIPseq_targets, `signalValue` ~ `symbol`, FUN = max, na.rm = TRUE) # max signalValue for gene
no_dup_ChIP <- rename(no_dup_ChIP, c('symbol' = 'SYMBOL'))
DEseq2_24h <- full_join(DEseq2_24h, no_dup_ChIP)

# merge RNA and prot
comb_prot_RNA <- merge(DEqMS_result, DEseq2_24h, by = "Protein.Ids")
comb_prot_RNA <- full_join(comb_prot_RNA, no_dup_ChIP)
comb_prot_RNA <- comb_prot_RNA %>% distinct(`gene`, .keep_all = TRUE)

RNA_diff <- comb_prot_RNA %>% dplyr::filter(abs(`logFC.y`) > 0.585 & `P-adj` < 0.05)
RNA_diff <- RNA_diff %>% distinct(`gene`, .keep_all = TRUE)

STAT1_prot <- comb_prot_RNA %>% dplyr::filter(`signalValue` > 0)
STAT1_prot <- STAT1_prot %>% distinct(`gene`, .keep_all = TRUE)

prot_diff <- comb_prot_RNA %>% dplyr::filter(abs(`logFC.x`) > 0.585 & `sca.adj.pval` < 0.05)
prot_diff <- prot_diff %>% distinct(`gene`, .keep_all = TRUE)

comb_prot_RNA$class <- rep("all proteins")
prot_diff$class <- rep("diff. nasc. prot.")
RNA_diff$class <- rep("diff. RNA abundance")
STAT1_prot$class <- rep("STAT1 targets")

comb_prot_RNA2 <- rbind(comb_prot_RNA, RNA_diff, STAT1_prot, prot_diff)
count_prots <- as.data.frame(cbind(`count` = c(nrow(comb_prot_RNA), nrow(RNA_diff), nrow(prot_diff), nrow(STAT1_prot)), 
                                   `class` = c("all proteins","diff. RNA abundance","diff. nasc. prot.","STAT1 targets")))

# plot scatter panel  
corr_panel <- ggplot() + 
  geom_abline(alpha = 0.6) + geom_hline(alpha = 0.6, yintercept = 0) + geom_vline(alpha = 0.6, xintercept = 0) +
  geom_point(data = comb_prot_RNA2, aes(x = `logFC.x`, y =  `logFC.y`, color = `class`), size = 3, show.legend = FALSE, alpha = 0.4) + facet_wrap(~class, ncol = 4) +
  stat_cor(data = comb_prot_RNA2, aes(label = ..r.label.., x = logFC.x, y =  `logFC.y`), r.accuracy = 0.001, method = "pearson", color = "black", size = 5, label.x = 1.4, label.y = -2) + 
  geom_text(data = count_prots, aes(x = 2.5, y =  5, label = paste0("n=",`count`)), size = 5) +
  xlab("log2(IFN-g/control) - nasc. prot.") + ylab("log2(IFN-g/control) - RNA-seq.") + theme_bw() + 
  theme(plot.title = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(color="black", size=14, face="bold"), axis.title.y = element_text(color="black", size=14, face="bold"), axis.text = element_text(size = 14),
        legend.title = element_text(color="black", size=14, face="bold"), legend.position = "bottom", legend.text = element_text(size = 14))

corr_panel


