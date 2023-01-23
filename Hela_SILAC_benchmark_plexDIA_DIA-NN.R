##Load the necessary packages
library(tidyverse)
library(BiocManager)
#library(RColorBrewer)
#library(ggalt)
#library(ggExtra)
#library(ggpubr)
library(cowplot)
#library(pheatmap)
#library(dendextend)
#library(ggthemes)
#library(Hmisc)
#library(PTXQC)
#library(Biostrings)
#library(matrixStats)
library(data.table)
#library(diann)
library(readxl)
library(iq)

#### set working directory !!!!!!!
setwd("~ your working directory and path to relevant files/")
getwd()

# import data and process
p_translated <- read_tsv(file = "report.pr_matrix_channels_ms1_translated.tsv") # MS1-based quantification
p_translated <- p_translated[!grepl("Cont_", p_translated$Protein.Ids),] # remove contaminants

# subset and label for each SILAC channel
p_translated_L <- p_translated %>% dplyr::select(`Protein.Group`:`Precursor.Id`, contains("raw-L"))
colnames(p_translated_L) = c("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description","Proteotypic","Stripped.Sequence","Modified.Sequence","Precursor.Charge","Precursor.Id",
                             paste0("AAAR1_DIA_L",seq(1:3)), paste0("AAAR2_DIA_L",seq(1:3)), paste0("AAAR1_MS1DIA_L",seq(1:3)), paste0("AAAR2_MS1DIA_L",seq(1:3)))

p_translated_M <- p_translated %>% dplyr::select(`Protein.Group`:`Precursor.Id`, contains("raw-M"))
colnames(p_translated_M) = c("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description","Proteotypic","Stripped.Sequence","Modified.Sequence","Precursor.Charge","Precursor.Id",
                             paste0("AAAR1_DIA_M",seq(1:3)), paste0("AAAR2_DIA_M",seq(1:3)), paste0("AAAR1_MS1DIA_M",seq(1:3)), paste0("AAAR2_MS1DIA_M",seq(1:3)))

p_translated_H <- p_translated %>% dplyr::select(`Protein.Group`:`Precursor.Id`, contains("raw-H"))
colnames(p_translated_H) = c("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description","Proteotypic","Stripped.Sequence","Modified.Sequence","Precursor.Charge","Precursor.Id",
                             paste0("AAAR1_DIA_H",seq(1:3)), paste0("AAAR2_DIA_H",seq(1:3)), paste0("AAAR1_MS1DIA_H",seq(1:3)), paste0("AAAR2_MS1DIA_H",seq(1:3)))


p_translated <- left_join(p_translated_L, p_translated_M) %>% left_join(p_translated_H)

pg_l <- p_translated %>% gather(contains(c("_L","_M","_H")), key = "Experiment", value = "precursor_translated")
pg_l$precursor_translated[pg_l$precursor_translated == 0] <- NA

SILAC_precursors <- pg_l[!is.na(pg_l$`precursor_translated`),] # remove 0 values for norm SILAC ratios


# add specified label info
prec_S <- SILAC_precursors %>% mutate(label = case_when(endsWith(Experiment, "_L1") ~ "L", endsWith(Experiment, "_L2") ~ "L", endsWith(Experiment, "_L3") ~ "L",
                                                        endsWith(Experiment, "_M1") ~ "M", endsWith(Experiment, "_M2") ~ "M", endsWith(Experiment, "_M3") ~ "M",
                                                        endsWith(Experiment, "_H1") ~ "H", endsWith(Experiment, "_H2") ~ "H", endsWith(Experiment, "_H3") ~ "H"))

# save prepared precursor tables
write_tsv(prec_S[,c("Protein.Ids","Genes","Precursor.Id","Experiment","precursor_translated")], "all_precursors_ldf_MS1.tsv")

# use iq R package to calculate protein level MaxLFQ values for the different SILAC channels
# MS1 quant - works better without median normalization in this benchmark
process_long_format("all_precursors_ldf_MS1.tsv", output_filename = "iq_pg_proc_MS1.tsv", primary_id = "Protein.Ids", secondary_id = c("Precursor.Id"),
                    annotation_col = c("Genes"), sample_id = "Experiment", intensity_col = "precursor_translated", normalization = "none", filter_double_less = NULL)

# import protein level quant data
prot_lfq <- read_tsv("iq_pg_proc_MS1.tsv")

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
prot_SF$`sum_int` <- log2(rowSums(2^(prot_SF[,c("int_H","int_M","int_L")]), na.rm = TRUE))
prot_SF$sum_int[prot_SF$sum_int == -Inf] <- NA

prot_SF <- prot_SF %>% mutate(perc = case_when(startsWith(Experiment, "AAAR1_DIA_") ~ "DIA m1", startsWith(Experiment, "AAAR2_DIA") ~ "DIA m1",
                                               startsWith(Experiment, "AAAR1_RS1DIA_") ~ "DIA m2", startsWith(Experiment, "AAAR2_RS1DIA") ~ "DIA m2")) # 

prot_SF <- prot_SF %>% mutate(sample = case_when(startsWith(Experiment, "AAAR1_") ~ "AAAR1", startsWith(Experiment, "AAAR2_") ~ "AAAR2")) # 


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

prot_int_agg <- aggregate(data = prot_SF, `sum_int` ~ Experiment + Protein.Ids  + `Genes` + perc + sample, FUN = median, na.rm = TRUE)

# !!! subset for method and SILAC ratio !!!
selected_method <- "DIA m1"

prot_HM_agg_f <- prot_HM_agg %>% dplyr::filter(`perc` == selected_method)
prot_HL_agg_f <- prot_HL_agg %>% dplyr::filter(`perc` == selected_method)
prot_ML_agg_f <- prot_ML_agg %>% dplyr::filter(`perc` == selected_method)
prot_int_agg_f <- prot_int_agg %>% dplyr::filter(`perc` == selected_method)

# create tables and save in complete collection! with different methods and quantification methods
prot_SILAC_DIA_m1_MS1 <- rbind(prot_HM_agg_f, prot_HL_agg_f, prot_ML_agg_f) %>% left_join(prot_int_agg_f)
prot_SILAC_DIA_m1_MS1$quant <- rep("MS1")

# !!! subset for method and SILAC ratio !!!
selected_method <- "DIA m2"

prot_HM_agg_f <- prot_HM_agg %>% dplyr::filter(`perc` == selected_method)
prot_HL_agg_f <- prot_HL_agg %>% dplyr::filter(`perc` == selected_method)
prot_ML_agg_f <- prot_ML_agg %>% dplyr::filter(`perc` == selected_method)
prot_int_agg_f <- prot_int_agg %>% dplyr::filter(`perc` == selected_method)

prot_SILAC_DIA_m2_MS1 <- rbind(prot_HM_agg_f, prot_HL_agg_f, prot_ML_agg_f) %>% left_join(prot_int_agg_f)
prot_SILAC_DIA_m2_MS1$quant <- rep("MS1")


# repeat processing for MS2 based quant data
p_translated <- read_tsv(file = "report.pr_matrix_channels_translated.tsv") # MS2-based quantification
p_translated <- p_translated[!grepl("Cont_", p_translated$Protein.Ids),] # remove contaminants

# subset and label for each SILAC channel
p_translated_L <- p_translated %>% dplyr::select(`Protein.Group`:`Precursor.Id`, contains("raw-L"))
colnames(p_translated_L) = c("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description","Proteotypic","Stripped.Sequence","Modified.Sequence","Precursor.Charge","Precursor.Id",
                             paste0("AAAR1_DIA_L",seq(1:3)), paste0("AAAR2_DIA_L",seq(1:3)), paste0("AAAR1_MS1DIA_L",seq(1:3)), paste0("AAAR2_MS1DIA_L",seq(1:3)))

p_translated_M <- p_translated %>% dplyr::select(`Protein.Group`:`Precursor.Id`, contains("raw-M"))
colnames(p_translated_M) = c("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description","Proteotypic","Stripped.Sequence","Modified.Sequence","Precursor.Charge","Precursor.Id",
                             paste0("AAAR1_DIA_M",seq(1:3)), paste0("AAAR2_DIA_M",seq(1:3)), paste0("AAAR1_MS1DIA_M",seq(1:3)), paste0("AAAR2_MS1DIA_M",seq(1:3)))

p_translated_H <- p_translated %>% dplyr::select(`Protein.Group`:`Precursor.Id`, contains("raw-H"))
colnames(p_translated_H) = c("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description","Proteotypic","Stripped.Sequence","Modified.Sequence","Precursor.Charge","Precursor.Id",
                             paste0("AAAR1_DIA_H",seq(1:3)), paste0("AAAR2_DIA_H",seq(1:3)), paste0("AAAR1_MS1DIA_H",seq(1:3)), paste0("AAAR2_MS1DIA_H",seq(1:3)))


p_translated <- left_join(p_translated_L, p_translated_M) %>% left_join(p_translated_H)

pg_l <- p_translated %>% gather(contains(c("_L","_M","_H")), key = "Experiment", value = "precursor_translated")
pg_l$precursor_translated[pg_l$precursor_translated == 0] <- NA

SILAC_precursors <- pg_l[!is.na(pg_l$`precursor_translated`),] # remove 0 values for norm SILAC ratios


# add specified label info
prec_S <- SILAC_precursors %>% mutate(label = case_when(endsWith(Experiment, "_L1") ~ "L", endsWith(Experiment, "_L2") ~ "L", endsWith(Experiment, "_L3") ~ "L",
                                                        endsWith(Experiment, "_M1") ~ "M", endsWith(Experiment, "_M2") ~ "M", endsWith(Experiment, "_M3") ~ "M",
                                                        endsWith(Experiment, "_H1") ~ "H", endsWith(Experiment, "_H2") ~ "H", endsWith(Experiment, "_H3") ~ "H"))

# save prepared precursor tables
write_tsv(prec_S[,c("Protein.Ids","Genes","Precursor.Id","Experiment","precursor_translated")], "all_precursors_ldf_MS2.tsv")


# use iq R package to calculate protein level MaxLFQ values for the different SILAC channels
# MS2 quant - works better with median normalization in this benchmark
process_long_format("all_precursors_ldf_MS2.tsv", output_filename = "iq_pg_proc_MS2.tsv", primary_id = "Protein.Ids", secondary_id = c("Precursor.Id"), annotation_col = c("Genes"), 
                    intensity_col_sep = ";", sample_id = "Experiment", intensity_col = "precursor_translated", normalization = "median", filter_double_less = NULL)

# import protein level quant data
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
prot_SF$`sum_int` <- log2(rowSums(2^(prot_SF[,c("int_H","int_M","int_L")]), na.rm = TRUE))
prot_SF$sum_int[prot_SF$sum_int == -Inf] <- NA

prot_SF <- prot_SF %>% mutate(perc = case_when(startsWith(Experiment, "AAAR1_DIA_") ~ "DIA m1", startsWith(Experiment, "AAAR2_DIA") ~ "DIA m1",
                                               startsWith(Experiment, "AAAR1_RS1DIA_") ~ "DIA m2", startsWith(Experiment, "AAAR2_RS1DIA") ~ "DIA m2")) # 

prot_SF <- prot_SF %>% mutate(sample = case_when(startsWith(Experiment, "AAAR1_") ~ "AAAR1", startsWith(Experiment, "AAAR2_") ~ "AAAR2")) # 

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

prot_int_agg <- aggregate(data = prot_SF, `sum_int` ~ Experiment + Protein.Ids  + `Genes` + perc + sample, FUN = median, na.rm = TRUE)

# !!! subset for method and SILAC ratio !!!
selected_method <- "DIA m1"

prot_HM_agg_f <- prot_HM_agg %>% dplyr::filter(`perc` == selected_method)
prot_HL_agg_f <- prot_HL_agg %>% dplyr::filter(`perc` == selected_method)
prot_ML_agg_f <- prot_ML_agg %>% dplyr::filter(`perc` == selected_method)
prot_int_agg_f <- prot_int_agg %>% dplyr::filter(`perc` == selected_method)

# create tables and save in complete collection! with different methods and quantification methods
prot_SILAC_DIA_m1_MS2 <- rbind(prot_HM_agg_f, prot_HL_agg_f, prot_ML_agg_f) %>% left_join(prot_int_agg_f)
prot_SILAC_DIA_m1_MS2$quant <- rep("MS2")

# !!! subset for method and SILAC ratio !!!
selected_method <- "DIA m2"

prot_HM_agg_f <- prot_HM_agg %>% dplyr::filter(`perc` == selected_method)
prot_HL_agg_f <- prot_HL_agg %>% dplyr::filter(`perc` == selected_method)
prot_ML_agg_f <- prot_ML_agg %>% dplyr::filter(`perc` == selected_method)
prot_int_agg_f <- prot_int_agg %>% dplyr::filter(`perc` == selected_method)

prot_SILAC_DIA_m2_MS2 <- rbind(prot_HM_agg_f, prot_HL_agg_f, prot_ML_agg_f) %>% left_join(prot_int_agg_f)
prot_SILAC_DIA_m2_MS2$quant <- rep("MS2")


## DDA from Maxquant
# import DDA ratio counts and values
prot_SILAC_DDA <- read_tsv(file = "Hela_SILAC_DDA_quant_proteins.tsv")

# stack all quant data
prot_SILAC_agg <- rbind(prot_SILAC_DIA_m1_MS1, prot_SILAC_DIA_m1_MS2, prot_SILAC_DIA_m2_MS1, prot_SILAC_DIA_m2_MS2, prot_SILAC_DDA)

# save and import complete quant data
prot_SILAC_agg <- unique(prot_SILAC_agg)
write_tsv(prot_SILAC_agg, "aggregated_SILAC_ratios_all_methods.tsv")
prot_SILAC_agg <- read_tsv(file ="aggregated_SILAC_ratios_all_methods.tsv")
prot_SILAC_agg <- prot_SILAC_agg[!is.na(prot_SILAC_agg$num_ratio),] # remove 0 values for norm SILAC ratios

prot_SILAC_agg$n <- rep(1)
prot_count <- aggregate(data = prot_SILAC_agg, `n` ~ `Experiment` + `perc` + `sample` + `ratio` + `quant`, FUN = sum, na.rm = TRUE)

# select sample
prot_SILAC_agg_f <- prot_SILAC_agg %>% dplyr::filter(`sample` == "AAAR1")
prot_count_f <- prot_count %>% dplyr::filter(`sample` == "AAAR1")

prot_count_mean <- aggregate(data = prot_count_f, `n` ~ perc + sample + ratio + quant, FUN = mean, na.rm = TRUE)

# plot
sample1 <- ggplot() + 
  geom_col(prot_count_mean, mapping = aes(x = `ratio`, fill = `quant`, y = n), stat = "count", position = "dodge", show.legend = FALSE, alpha = 0.9) +
  geom_point(prot_count_f, mapping = aes(x = `ratio`, group = `quant`, shape = `quant`, y = `n`), position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), show.legend = FALSE, size = 3, alpha = 0.6) +
  geom_text(prot_count_mean, mapping = aes(x = `ratio`, group = `quant`, label = round(n, digits = 0), y = n), vjust = 0.5, size = 5, 
            color = "black", position = position_dodge(width = 0.95), angle = 90, hjust = 1.5) + 
  facet_wrap(~`perc`) + ylab("quantified proteins") + xlab("") + theme_bw() + ylim(0,7000) + ggtitle("SILAC mix 1") + scale_fill_manual(values = c("dodgerblue3","firebrick3")) +
  theme(plot.title = element_text(color="black", size = 16, face= "bold", vjust = 0.5), axis.title.x = element_text(color="black", size=14, face= "bold"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color="black", size = 14), axis.text.y = element_text(color="black", size=14), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(), axis.title.y = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 14, face = "bold"), legend.position = "right",
        legend.text = element_text(color="black", size=14), legend.title = element_text(color="black", size=14, face="bold"))

# select sample
prot_SILAC_agg_f <- prot_SILAC_agg %>% dplyr::filter(`sample` == "AAAR2")
prot_count_f <- prot_count %>% dplyr::filter(`sample` == "AAAR2")

prot_count_mean <- aggregate(data = prot_count_f, `n` ~ perc + sample + ratio + quant, FUN = mean, na.rm = TRUE)

# plot
sample2 <- ggplot() + 
  geom_col(prot_count_mean, mapping = aes(x = `ratio`, fill = `quant`, y = n), stat = "count", position = "dodge", show.legend = FALSE, alpha = 0.9) +
  geom_point(prot_count_f, mapping = aes(x = `ratio`, group = `quant`, shape = `quant`, y = `n`), position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), show.legend = FALSE, size = 3, alpha = 0.6) +
  geom_text(prot_count_mean, mapping = aes(x = `ratio`, group = `quant`, label = round(n, digits = 0), y = n), vjust = 0.5, size = 5, 
            color = "black", position = position_dodge(width = 0.95), angle = 90, hjust = 1.5) +
  facet_wrap(~`perc`) + ylab("quantified proteins") + xlab("") + theme_bw() + ylim(0,7000) + ggtitle("SILAC mix 2") + scale_fill_manual(values = c("dodgerblue3","firebrick3")) +
  theme(plot.title = element_text(color="black", size = 16, face= "bold", vjust = 0.5), axis.title.x = element_text(color="black", size=14, face= "bold"),
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color="black", size = 14), axis.text.y = element_blank(), panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(), axis.title.y = element_blank(), strip.text = element_text(size = 14, face = "bold"), legend.position = "right",
      legend.text = element_text(color="black", size=14), legend.title = element_text(color="black", size=14, face="bold"))

plot_grid(sample1, sample2, ncol = 2, rel_widths = c(0.55,0.45))
top <- plot_grid(sample1, sample2, ncol = 2, rel_widths = c(0.55,0.45))
 

# calculate additional quant metrics
sd_l <- aggregate(data = prot_SILAC_agg, `num_ratio` ~ `perc` + `Genes` + `Protein.Ids` + `sample` + `quant` + `ratio`, FUN = sd, na.rm = TRUE) # select corret ratio!
colnames(sd_l) <- c("perc","Gene_names","Protein.Group","sample","quant","ratio","sd")

mean_l <- aggregate(data = prot_SILAC_agg, `num_ratio` ~ `perc` + `Genes` + `Protein.Ids` + `sample` + `quant` + `ratio`, FUN = mean, na.rm = TRUE) # select corret ratio!
colnames(mean_l) <- c("perc","Gene_names","Protein.Group","sample","quant","ratio","mean")

stats_l <- merge(mean_l, sd_l)
stats_l$cv <- (stats_l$sd / stats_l$mean)*100

# save for respective SIKLAC ratio
write_tsv(stats_l, "proteinGroups_stats_l_complete.tsv")
stats_l <- read_tsv(file = "proteinGroups_stats_l_complete.tsv")

# CV 
stats_lf <- stats_l %>% dplyr::filter(`sample` == "AAAR1") # select indivdual sample for plot and plotgrid

p_meds <- aggregate(list(stats_lf$`cv`), by = list(stats_lf$`perc`, stats_lf$`sample`, stats_lf$`ratio`, stats_lf$`quant`), median, na.rm = TRUE)
colnames(p_meds) <- c("perc","sample","ratio","quant","cv")

AAAR1 <- ggplot(data = stats_lf, aes(y = `cv`, x = `ratio`, fill = `quant`)) +
  geom_boxplot(show.legend = FALSE, alpha = 0.9, outlier.alpha = 0.1) + #scale_y_continuous(limits = c(0,1)) +
  geom_text(data = p_meds, aes(y = -25, x = `ratio`, group = `quant`, label = round(`cv`, digits = 1)), col = "black", alpha = 0.9, size = 5, vjust = 0.5, 
            color = "black", position = position_dodge(width = 0.95), angle = 90) + 
  theme_bw() + facet_grid(~perc) + ggtitle("SILAC mix 1") +
  xlab("") + ylab("coefficient of variation [%]") + ylim(c(-45,150)) + scale_fill_manual(values = c("dodgerblue3","firebrick3")) +
  theme(plot.title = element_text(color="black", size = 14, face= "bold", vjust = 0.5), axis.title.x = element_text(color="black", size=14, face= "bold"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color="black", size = 14), axis.text.y = element_text(color="black", size=14), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(), axis.title.y = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 14, face = "bold"), legend.position = "right",
        legend.text = element_text(color="black", size=14), legend.title = element_text(color="black", size=14, face="bold"))


stats_lf <- stats_l %>% dplyr::filter(`sample` == "AAAR2") # select indivdual sample for plot and plotgrid

p_meds <- aggregate(list(stats_lf$`cv`), by = list(stats_lf$`perc`, stats_lf$`sample`, stats_lf$`ratio`, stats_lf$`quant`), median, na.rm = TRUE)
colnames(p_meds) <- c("perc","sample","ratio","quant","cv")

AAAR2 <- ggplot(data = stats_lf, aes(y = `cv`, x = `ratio`, fill = `quant`)) +
  geom_boxplot(show.legend = FALSE, alpha = 0.9, outlier.alpha = 0.1) + #scale_y_continuous(limits = c(0,1)) +
  geom_text(data = p_meds, aes(y = -25, x = `ratio`, group = `quant`, label = round(`cv`, digits = 1)), col = "black", alpha = 0.9, size = 5, vjust = 0.5, 
            color = "black", position = position_dodge(width = 0.95), angle = 90) + 
  theme_bw() + facet_grid(~perc) + xlab("") + ylab("") + 
  ylim(c(-45,150)) + scale_fill_manual(values = c("dodgerblue3","firebrick3")) + ggtitle("SILAC mix 2") +
  theme(plot.title = element_text(color="black", size = 14, face= "bold", vjust = 0.5), axis.title.y = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_text(color="black", size=14), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(), axis.title.x = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 14, face = "bold"), legend.position = "right",
        legend.text = element_text(color="black", size=14), legend.title = element_text(color="black", size=14, face="bold"))

plot_grid(AAAR1, AAAR2, ncol = 2, rel_widths = c(0.55,0.45))
middle <- plot_grid(AAAR1, AAAR2, ncol = 2, rel_widths = c(0.55,0.45))


### box plots of SILAC ratios to check quant accuracy
prot_SILAC_agg_f2 <- prot_SILAC_agg %>% dplyr::filter(`sample` == "AAAR1")
prot_SILAC_agg_f2 <- aggregate(data = prot_SILAC_agg_f2, `num_ratio` ~ `perc` + `Genes` + `Protein.Ids` + `sample` + `quant` + `ratio`, FUN = mean, na.rm = TRUE) # calculate mean value of SILAC ratio

# calculate difference from expected ratio
prot_SILAC_agg_f2 <- prot_SILAC_agg_f2 %>% mutate(expected = case_when(startsWith(ratio, "H/L") ~ (15/70), startsWith(ratio, "H/M") ~ 1, startsWith(ratio, "M/L") ~ (15/70))) # add numeric percentage to samples # step 2
prot_SILAC_agg_f2$log2error <- (log2(prot_SILAC_agg_f2$num_ratio) - log2(prot_SILAC_agg_f2$expected))

mix1_error <- aggregate(data = prot_SILAC_agg_f2, abs(`log2error`) ~ `perc` + `sample` + `quant` + `ratio`, FUN = median, na.rm = TRUE) # select corret ratio!

mix1 <- ggplot(data = prot_SILAC_agg_f2, aes(y = log2(`num_ratio`), x = `ratio`, fill = `quant`)) +
  geom_hline(yintercept = c(-2.222392,0), alpha = 0.5, size = 0.5, color = "black") +
  geom_boxplot(show.legend = FALSE, alpha = 0.9, outlier.alpha = 0.1) + scale_y_continuous(limits = c(-6,4)) +
  geom_text(data = mix1_error, aes(y = -5.5, x = `ratio`, group = `quant`, label = round(`abs(log2error)`, digits = 2)), col = "black", alpha = 0.9, size = 5, vjust = 0.5, 
            color = "black", position = position_dodge(width = 0.95), angle = 90) + 
  theme_bw() + facet_grid(~perc) + ggtitle("SILAC mix 1") + xlab("") + ylab("log2(SILAC ratio)") + scale_fill_manual(values = c("dodgerblue3","firebrick3")) +
  theme(plot.title = element_text(color="black", size = 14, face= "bold", vjust = 0.5), axis.title.x = element_text(color="black", size=14, face= "bold"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color="black", size = 14), axis.text.y = element_text(color="black", size=14), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(), axis.title.y = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 14, face = "bold"), legend.position = "right",
        legend.text = element_text(color="black", size=14), legend.title = element_text(color="black", size=14, face="bold"))


prot_SILAC_agg_f2 <- prot_SILAC_agg %>% dplyr::filter(`sample` == "AAAR2")
prot_SILAC_agg_f2 <- aggregate(data = prot_SILAC_agg_f2, `num_ratio` ~ `perc` + `Genes` + `Protein.Ids` + `sample` + `quant` + `ratio`, FUN = mean, na.rm = TRUE) # calculate mean value of SILAC ratio

prot_SILAC_agg_f2 <- prot_SILAC_agg_f2 %>% mutate(expected = case_when(startsWith(ratio, "H/L") ~ 2, startsWith(ratio, "H/M") ~ 1, startsWith(ratio, "M/L") ~ 2)) # add numeric percentage to samples # step 2
prot_SILAC_agg_f2$log2error <- (log2(prot_SILAC_agg_f2$num_ratio) - log2(prot_SILAC_agg_f2$expected))

mix2_error <- aggregate(data = prot_SILAC_agg_f2, abs(`log2error`) ~ `perc` + `sample` + `quant` + `ratio`, FUN = median, na.rm = TRUE) # select corret ratio!

mix2 <- ggplot(data = prot_SILAC_agg_f2, aes(y = log2(`num_ratio`), x = `ratio`, fill = `quant`)) +
  geom_hline(yintercept = c(1,0), alpha = 0.5, size = 0.5, color = "black") +
  geom_boxplot(show.legend = FALSE, alpha = 0.9, outlier.alpha = 0.1) + scale_y_continuous(limits = c(-6,4)) +
  geom_text(data = mix2_error, aes(y = -5.5, x = `ratio`, group = `quant`, label = round(`abs(log2error)`, digits = 2)), col = "black", alpha = 0.9, size = 5, vjust = 0.5, 
            color = "black", position = position_dodge(width = 0.95), angle = 90) + 
  theme_bw() + facet_grid(~perc) + xlab("") + ylab("") + scale_fill_manual(values = c("dodgerblue3","firebrick3")) + ggtitle("SILAC mix 2") +
  theme(plot.title = element_text(color="black", size = 14, face= "bold", vjust = 0.5), axis.title.y = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_text(color="black", size=14), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(), axis.title.x = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 14, face = "bold"), legend.position = "right",
        legend.text = element_text(color="black", size=14), legend.title = element_text(color="black", size=14, face="bold"))

# plot
plot_grid(mix1, mix2, ncol = 2, rel_widths = c(0.55,0.45))
bottom <- plot_grid(mix1, mix2, ncol = 2, rel_widths = c(0.55,0.45))

### plot method schemes
methodscheme <- read_xlsx("method_schemes_for_visualization.xlsx")
cycle_time <- read_xlsx("cycle_time.xlsx")

meth <- ggplot(data = methodscheme) + 
  geom_rect(aes(fill = `type`, xmin = startmz, xmax = endmz, ymin = iw-0.3, ymax = iw+0.3), alpha = 0.9, size = 8) + facet_wrap(~method, ncol = 3, scales = "free_x") +
  geom_text(data = cycle_time, mapping = aes(x = mz, y = iw, label = `cycle_time`), size = 5) + 
  xlab("m/z") + ylab("scan number") + theme_bw() + scale_x_continuous(breaks = c(300,600,900,1200)) + scale_fill_manual(values = c("dodgerblue3","firebrick3")) +
  theme(plot.title = element_text(color="black", size=12, face= "bold"), axis.title.x = element_text(color="black", size=14, face= "bold"),
        axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(color="black", size=14), 
        axis.title.y = element_text(color="black", size=14, face="bold"), legend.position = c(0.45,0.80), strip.text = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14), legend.title = element_blank(), legend.background = element_blank())


samplecomp <- read_xlsx("AAAR12_sample_composition.xlsx")

samp <- ggplot(data = samplecomp, aes(fill = `channel`, x = `sample`, y = `composition`, label = `label`), ) + 
  geom_col(alpha = 0.9, show.legend = FALSE) + 
  geom_text(position = position_stack(vjust = .5), size = 5) + scale_fill_brewer(palette = "Set2") +
  xlab("") + ylab("labeled Hela lysate [µg]") + theme_bw() + theme(plot.title = element_text(color="black", size=12, face= "bold"), axis.title.x = element_text(color="black", size=14, face= "bold"),
                                                                   axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(color="black", size=14), 
                                                                   axis.title.y = element_text(color="black", size=14, face="bold"), legend.position = c(0.45,0.80), strip.text = element_text(size = 14, face = "bold"),
                                                                   legend.text = element_text(size = 14), legend.title = element_blank(), legend.background = element_blank())


# plot final full panel (figure 4 in manuscript)
SampMeth <- plot_grid(samp, meth, ncol = 2, rel_widths = c(0.4,0.6))
plot_grid(SampMeth, top, middle, bottom, ncol = 1, rel_heights = c(0.2,0.2,0.2,0.3))

