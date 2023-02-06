##Load the necessary packages
library(tidyverse)
library(BiocManager)
library(ggrepel)
library(reshape2)
library(naniar)
library(RColorBrewer)
library(ggalt)
library(ggExtra)
library(ggpubr)
library(cowplot)
library(Biostrings)
library(Biostrings)
library(readxl)
library(Peptides)

#### set working directory !!!!!!!
setwd("~ path to your working directory")
getwd()

# color palette for plots
nb.cols <- 9
greens <- colorRampPalette(brewer.pal(8, "Greens"))(nb.cols)
greens <- c(greens,"#E08B00")
greens <- greens[-c(8,9)]
  
nb.cols <- 15
blues <- colorRampPalette(brewer.pal(8, "Blues"))(nb.cols)

## import txt file output tables from Maxquant
# import evidence table (can be obtained from PRIDE repository with ID of dataset PXD036886)
evidence <- read_tsv(file = "evidence.txt", na = "NaN", col_types = cols(`Potential contaminant` = col_character(), `Reverse` = col_character()))
colnames(evidence) <- str_replace_all(colnames(evidence), "\\s",replacement = "_")
colnames(evidence)

evidence <- evidence %>% dplyr::filter(`Potential_contaminant` != "+" | `Reverse` != "+") 
evidence <- evidence[!grepl("CON_", evidence$Proteins),] # remove contaminants
evidence <- evidence[!grepl("REV_", evidence$Leading_proteins),] # remove decoys

evidence <- evidence %>% mutate(input = case_when(startsWith(Experiment, "1ug") ~ "1", startsWith(Experiment, "5ug") ~ "5", startsWith(Experiment, "10ug") ~ "10",
                                                      startsWith(Experiment, "25ug") ~ "25", startsWith(Experiment, "50ug") ~ "50", startsWith(Experiment, "100ug") ~ "100",
                                                      startsWith(Experiment, "150ug") ~ "150", startsWith(Experiment, "200ug") ~ "200", startsWith(Experiment, "250ug") ~ "250",
                                                      startsWith(Experiment, "300ug") ~ "300", startsWith(Experiment, "notEnriched") ~ "nE",
                                                      startsWith(Experiment, "2uL") ~ "2", startsWith(Experiment, "3uL") ~ "3",
                                                      startsWith(Experiment, "4uL") ~ "4", startsWith(Experiment, "5uL") ~ "5", startsWith(Experiment, "6uL") ~ "6",
                                                      startsWith(Experiment, "7uL") ~ "7", startsWith(Experiment, "8uL") ~ "8"))

# import and process protein groups table (can be obtained from PRIDE repository with ID of dataset PXD036886)
prot <- read_tsv(file = "proteinGroups.txt", na = "NaN", col_types = cols(Reverse = col_character(), `Only identified by site` = col_character()))

# remove spaces in column names 
colnames(prot) <- str_replace_all(colnames(prot), "\\s",replacement = "_")
colnames(prot)

# Remove contaminants, reverse hits and only identified by site
prot_f <- prot %>% filter(Only_identified_by_site != "+", Reverse != "+", Potential_contaminant != "+")
prot_f <- prot_f[!grepl("CON_", prot_f$Protein_IDs),] # remove contaminants
prot_f <- prot_f[!grepl("REV_", prot_f$Protein_IDs),] # remove decoys

# Select columns that we will need for further processing
prot_f1 <- prot_f %>% dplyr::select(contains("Protein"), Protein_IDs:Number_of_proteins,
                                    starts_with("Peptides_"), matches("^Sequence_coverage_[^[]"), `Mol._weight_[kDa]`,
                                    starts_with("Identification"), matches("Ratio_./._[^vit]"), matches("iBAQ_"), matches("^Intensity_._."))

# create smaller long data frames with focus on various parameters of protein groups
# peptides table
Peptides_tb <- prot_f1 %>% dplyr::select(Protein_IDs:Peptides_notEnriched_R2) %>% 
  gather(Peptides_100ug_R1:Peptides_notEnriched_R2 ,key = "Experiment", value = "Peptide_Number") %>%
  mutate(Experiment = str_remove_all(Experiment,pattern = "Peptides_"))
Peptides_tb[(Peptides_tb) == 0] <- NA

# sequence coverage
Seq_cov_tb <- prot_f1 %>% dplyr::select(Protein_IDs:Fasta_headers,starts_with("Seq")) %>% 
  gather(starts_with("Seq"), key = "Experiment", value = "Seq_cov_[%]") %>%
  mutate(Experiment = str_remove_all(Experiment,pattern = "Sequence_coverage_")) %>%
  mutate(Experiment = str_remove_all(Experiment,pattern ="_\\[%]"))

# SILAC-ratio count
ratio_count_tb <- prot_f1 %>% dplyr::select(Protein_IDs:Fasta_headers, starts_with("Ratio_H/M_count")) %>% 
  gather(starts_with("Ratio_H/M_count_"), key = "Experiment", value = "Ratio_H/M_count") %>%
  mutate(Experiment = str_remove_all(Experiment, pattern = "Ratio_H/M_count_"))
ratio_count_tb[(ratio_count_tb) == 0] <- NA

# intensities
Intensity_L <- prot_f1 %>% dplyr::select(Protein_IDs:Fasta_headers,starts_with("Intensity")) %>% 
  gather(starts_with("Intensity_L"), key = "Experiment", value = "Intensity_L") %>% 
  mutate(Experiment = str_remove_all(Experiment, pattern = "Intensity_L_")) %>%
  dplyr::select(Protein_IDs:Fasta_headers, Experiment:Intensity_L)

Intensity_M <-  prot_f1 %>% dplyr::select(Protein_IDs:Fasta_headers,starts_with("Intensity")) %>% 
  gather(starts_with("Intensity_M"), key = "Experiment", value = "Intensity_M") %>% 
  mutate(Experiment = str_remove_all(Experiment, pattern = "Intensity_M_")) %>% 
  dplyr::select(Protein_IDs:Fasta_headers, Experiment:Intensity_M)

Intensity_H <- prot_f1 %>% dplyr::select(Protein_IDs:Fasta_headers,starts_with("Intensity")) %>% 
  gather(starts_with("Intensity_H"), key = "Experiment", value = "Intensity_H") %>% 
  mutate(Experiment = str_remove_all(Experiment, pattern = "Intensity_H_")) %>% 
  dplyr::select(Protein_IDs:Fasta_headers, Experiment:Intensity_H)

# combine intensities
Intensity_table <- left_join(Intensity_M,Intensity_H) %>% left_join(Intensity_L)
Intensity_table[(Intensity_table) == 0] <- NA

# normalized SILAC ratios
# H/L
Silac_tb_norm_HL <- prot_f1 %>% dplyr::select(Protein_IDs:Fasta_headers,starts_with("Ratio"),-("Ratio_M/L_normalized":"Ratio_H/M_count")) %>% 
  gather(contains("H/L_normalized"), key = "Experiment", value = "Ratio_norm_H/L") %>% 
  mutate(Experiment = str_remove_all(Experiment, pattern = "Ratio_H/L_normalized_")) %>%
  dplyr::select(Protein_IDs:Fasta_headers, Experiment:"Ratio_norm_H/L")

# M/L
Silac_tb_norm_ML <- prot_f1 %>%
  dplyr::select(Protein_IDs:Fasta_headers,starts_with("Ratio"),-("Ratio_M/L_normalized":"Ratio_H/M_count")) %>% 
  gather(contains("M/L_normalized"), key = "Experiment", value = "Ratio_norm_M/L") %>% 
  mutate(Experiment = str_remove_all(Experiment,pattern = "Ratio_M/L_normalized_")) %>%
  dplyr::select(Protein_IDs:Fasta_headers,Experiment:"Ratio_norm_M/L")

# H/M
Silac_tb_norm_HM <- prot_f1 %>% dplyr::select(Protein_IDs:Fasta_headers,starts_with("Ratio"), -("Ratio_M/L_normalized":"Ratio_H/M_count")) %>% 
  gather(contains("H/M_normalized"), key = "Experiment", value = "Ratio_norm_H/M") %>% 
  mutate(Experiment = str_remove_all(Experiment, pattern = "Ratio_H/M_normalized_")) %>%
  dplyr::select(Protein_IDs:Fasta_headers,Experiment:"Ratio_norm_H/M")

# combine normalized SILAC ratios
Silac_tb_norm <- left_join(Silac_tb_norm_HL,Silac_tb_norm_ML) %>% left_join(Silac_tb_norm_HM)

# unnormalized ratios
Silac_tb_unnorm_HM <- prot_f1 %>% dplyr::select(Protein_IDs:Fasta_headers, starts_with("Ratio_H/M_"), -starts_with(c("Ratio_H/M_normalized","Ratio_H/M_normalized_","Ratio_H/M_count_","Ratio_H/M_count"))) %>% 
  gather(contains("Ratio_H/M_"), key = "Experiment", value = "Ratio_H/M") %>% 
  mutate(Experiment = str_remove_all(Experiment, pattern = "Ratio_H/M_")) %>%
  dplyr::select(Protein_IDs:Fasta_headers, Experiment:`Ratio_H/M`)

Silac_tb_unnorm_HL <- prot_f1 %>% dplyr::select(Protein_IDs:Fasta_headers, starts_with("Ratio_H/L_"), -starts_with(c("Ratio_H/L_normalized","Ratio_H/L_normalized_","Ratio_H/L_count_","Ratio_H/L_count"))) %>% 
  gather(contains("Ratio_H/L_"), key = "Experiment", value = "Ratio_H/L") %>% 
  mutate(Experiment = str_remove_all(Experiment, pattern = "Ratio_H/L_")) %>%
  dplyr::select(Protein_IDs:Fasta_headers, Experiment:`Ratio_H/L`)

Silac_tb_unnorm_ML <- prot_f1 %>% dplyr::select(Protein_IDs:Fasta_headers, starts_with("Ratio_M/L_"), -starts_with(c("Ratio_M/L_normalized","Ratio_M/L_normalized_","Ratio_M/L_count_","Ratio_M/L_count"))) %>% 
  gather(contains("Ratio_M/L_"), key = "Experiment", value = "Ratio_M/L") %>% 
  mutate(Experiment = str_remove_all(Experiment, pattern = "Ratio_M/L_")) %>%
  dplyr::select(Protein_IDs:Fasta_headers, Experiment:`Ratio_M/L`)

# Create combined long table
table_merg <- left_join(Silac_tb_norm, Silac_tb_unnorm_HM) %>% left_join(Silac_tb_unnorm_HL) %>% left_join(Silac_tb_unnorm_ML) %>% left_join(Intensity_table) %>% left_join(Peptides_tb) %>% left_join(Seq_cov_tb)

table_merg_f <- table_merg
table_merg_f <- table_merg_f %>% mutate(input = case_when(startsWith(Experiment, "1ug") ~ "1", startsWith(Experiment, "5ug") ~ "5", startsWith(Experiment, "10ug") ~ "10",
                                                          startsWith(Experiment, "25ug") ~ "25", startsWith(Experiment, "50ug") ~ "50", startsWith(Experiment, "100ug") ~ "100",
                                                          startsWith(Experiment, "150ug") ~ "150", startsWith(Experiment, "200ug") ~ "200", startsWith(Experiment, "250ug") ~ "250",
                                                          startsWith(Experiment, "300ug") ~ "300", startsWith(Experiment, "notEnriched") ~ "nE",
                                                          startsWith(Experiment, "2uL") ~ "2", startsWith(Experiment, "3uL") ~ "3",
                                                          startsWith(Experiment, "4uL") ~ "4", startsWith(Experiment, "5uL") ~ "5", startsWith(Experiment, "6uL") ~ "6",
                                                          startsWith(Experiment, "7uL") ~ "7", startsWith(Experiment, "8uL") ~ "8"))

table_merg_f$SILAC_enrich <- (table_merg_f$`Ratio_H/L` + table_merg_f$`Ratio_M/L`)

# save and import table
write_tsv(table_merg_f, "table_merg_f.tsv", na = "NA")
table_merg_f <- read_tsv(file = "table_merg_f.tsv")



## continue with analysis of the 2 different data sets separately ##

#### analysis of data with MAA bead dilution series and contrast to sample without enrichment 
# subset evidence precursor level data to exclude protein input dilution data in this section
evidence_f <- evidence[!grepl("ug_", evidence$Experiment),] # remove input dilution samples
colnames(evidence_f)

ggplot(evidence_f, aes(x = `Retention_time`, fill = `Experiment`)) +
  geom_histogram(show.legend = FALSE, linewidth = 1, alpha = 0.2, bins = 105, color = "black") + #geom_smooth(show.legend = FALSE, color = "black", size = 1.5, alpha = 0.7) +
  scale_y_continuous(limits = c()) + facet_wrap('Experiment', ncol = 2) +
  theme_bw() + xlab("Retention time [min]") + ylab("identified precursors per min") +
  theme(axis.title.x = element_text(size = 14, face = "bold"), axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14, face = "bold"), axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold"))

# calculate S/N ratio
evidence_f$SILAC_enrich <- (evidence_f$`Ratio_H/L` + evidence_f$`Ratio_M/L`)
evidence_f$n <- rep(1)

# MAA dilution
neworder <- c("2","3","4","5","6","7","8","nE")
evidence_f <- evidence_f %>% arrange(transform(evidence_f, `input` = factor(`input`, levels = neworder)), `input`)

# plot SILAC enrich ratios
p_meds <- aggregate(evidence_f$`SILAC_enrich`, by = list(evidence_f$`input`), FUN = median, na.rm = TRUE)

SILAC_violin <- ggplot(data = evidence_f, aes(x = as.factor(`input`), y = log2(`SILAC_enrich`))) +
  geom_hline(yintercept = 0, linetype="dashed", color = "black", size = 0.25) +
  geom_violin(aes(fill = as.factor(`input`)), show.legend = FALSE) + geom_boxplot(aes(fill = as.factor(`input`)), width = 0.2, show.legend = FALSE, outlier.size = -1) + 
  theme_bw() + xlab("bead volume [無]") + ylab("log2(labeled-/unlabeled\nprecursor intensity)") + scale_y_continuous(limits = c(min(log2(evidence_f$SILAC_enrich)), max(log2(evidence_f$SILAC_enrich), na.rm = T)+2)) +
  scale_fill_manual(values = greens) +
  geom_text(data = p_meds, aes(x = as.factor(`Group.1`), y = max(log2(evidence_f$SILAC_enrich), na.rm = T)+1, label = round(log2(`x`), digits = 2)), size = 5, vjust = 0.5, angle = 0) +
  theme(plot.title = element_text(color="black", size=12, face= "bold"), axis.title.x = element_text(color="black", size=14, face= "bold"),
        axis.text.x = element_text(color = "black", size=14, angle = 0),
        axis.text.y = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14, face="bold"))
SILAC_violin

prec_HM <- evidence_f[,c("Experiment","input","n","Ratio_H/M")]
prec_HM[(prec_HM) == 0] <- NA
prec_HM <- prec_HM[!is.na(prec_HM$`Ratio_H/M`),] # remove NA values

all_precursors_agg <- aggregate(data = prec_HM, `n` ~ `Experiment` + `input`, FUN = sum, na.rm = TRUE)
all_precursors_agg$class <- rep("quantified")
all_precursors_agg$class <- rep("a")

p_meds <- aggregate(list(all_precursors_agg$`n`), by = list(all_precursors_agg$`input`), median)
colnames(p_meds) <- c("input","n")
p_meds$class <- rep("quantified")

quant_prec <- ggplot() +
  geom_col(data = p_meds, aes(x = as.factor(`input`), y = `n`, fill = as.factor(`input`)), color = "black", show.legend = FALSE) +
  geom_jitter(data = all_precursors_agg, aes(x = as.factor(`input`), y = `n`), position = position_jitter(width = 0.25), show.legend = FALSE, size = 4, alpha = 0.4) +
  geom_text(data = p_meds, aes(x = as.factor(`input`), label = round(n, digits = 0), y = n), vjust = 0.5, hjust = 2.0, size = 5, color = "black", angle = 90) +
  ylab("quantified precursors") + xlab("") + theme_bw() +  scale_fill_manual(values = greens) +
  theme(plot.title = element_text(color="black", size=16, face= "bold"), axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(vjust = 1, color = "black", size=14),
        axis.text.y = element_text(color = "black", size = 14), axis.title.y = element_text(color = "black", size = 14, face = "bold"))
quant_prec

# identified _precursors
prec_int <- evidence_f[,c("Experiment","input","n","Intensity")]
prec_int[(prec_int) == 0] <- NA
prec_int <- prec_int[!is.na(prec_int$Intensity),] # remove NA values

all_ident_precursors_agg <- aggregate(data = prec_int, `n` ~ `Experiment` + `input`, FUN = sum, na.rm = TRUE)
all_ident_precursors_agg$class <- rep("b")

p_meds_ident <- aggregate(list(all_ident_precursors_agg$`n`), by = list(all_ident_precursors_agg$`input`, all_ident_precursors_agg$`class`), median)
colnames(p_meds_ident) <- c("input","class","n")

p_meds <- aggregate(list(all_precursors_agg$`n`), by = list(all_precursors_agg$`input`), median)
colnames(p_meds) <- c("input","n")
p_meds$class <- rep("a")

all_precursors_agg_neg <- all_precursors_agg
all_precursors_agg_neg$n <- all_precursors_agg_neg$n * -1

all_ident_precursors_diff <- rbind(all_ident_precursors_agg, all_precursors_agg_neg)
all_ident_precursors_diff <- aggregate(data = all_ident_precursors_diff, `n` ~ `Experiment` + `input`, FUN = sum, na.rm = TRUE)
all_ident_precursors_diff$class <- rep("b")

stack_val <- rbind(all_ident_precursors_diff, all_precursors_agg)
stack_val <- aggregate(data = stack_val, `n` ~ `input` + `class`, FUN = median, na.rm = TRUE)

stack_val2 <- rbind(all_precursors_agg, all_ident_precursors_agg)
stack_val2 <- aggregate(data = stack_val2, `n` ~ `input` + `class`, FUN = median, na.rm = TRUE)

full_meds <- rbind(p_meds, p_meds_ident)

jitter_full <- rbind(all_precursors_agg, all_ident_precursors_agg)

# plot identified (on top) and quantified precurosors (on bottom)
ident_prec <- ggplot() +
  geom_col(data = stack_val, aes(x = as.factor(`input`), y = `n`, fill = as.factor(`input`)), color = "black", show.legend = FALSE) +
  geom_jitter(data = jitter_full, aes(x = as.factor(`input`), y = `n`, shape = `class`), position = position_jitter(width = 0.25), show.legend = FALSE, size = 3, alpha = 0.4) +
  geom_text(data = full_meds %>% dplyr::filter(`class` == "a"), aes(x = as.factor(`input`), label = round(n, digits = 0), y = 10000), vjust = 0.5, size = 5, color = "black", angle = 90) +
  geom_text(data = full_meds %>% dplyr::filter(`class` == "b"), aes(x = as.factor(`input`), label = round(n, digits = 0), y = 25750), vjust = 0.5, size = 5, color = "black", angle = 90) +
  ylab("number of precursors") + xlab("bead volume [無]") + theme_bw() +
  scale_fill_manual(values = greens) +
  theme(plot.title = element_text(color="black", size=16, face= "bold"), axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(vjust = 0.5, color = "black", size=14, angle = 0, hjust = 0.5),
        axis.text.y = element_text(color = "black", size = 14), axis.title.y = element_text(color = "black", size = 14, face = "bold"))
ident_prec

## subset protein group, protein level data to exclude protein input dilution data in this section
table_merg_f <- table_merg_f[!grepl("ug_", table_merg_f$Experiment),] # remove input dilution samples

# number of quantified proteins per sample
table_merg_f2 <- table_merg_f[!is.na(table_merg_f$`Ratio_norm_H/M`),] # remove 0 values for norm SILAC ratios

# 
prot_n <- aggregate(table_merg_f2$`Ratio_norm_H/M`, by = list(table_merg_f2$`Experiment`,table_merg_f2$`input`), FUN = length)
colnames(prot_n) <- c("Experiment","input","quantified_proteins")

# plot number of quantified proteins in bar chart
table_merg_f2$n <- rep(1)
all_prot_agg <- aggregate(data = table_merg_f2, `n` ~ `Experiment` + `input`, FUN = sum, na.rm = TRUE)

p_meds <- aggregate(list(all_prot_agg$`n`), by = list(all_prot_agg$`input`), median)
colnames(p_meds) <- c("input","n")
p_meds$group <- rep("quant")

quant_prot <- ggplot() +
  geom_col(data = p_meds, aes(x = as.factor(`input`), y = n, fill = as.factor(`input`)), color = "black", position = "dodge", show.legend = FALSE) +
  geom_jitter(data = prot_n, aes(x = as.factor(`input`), y = `quantified_proteins`), position = position_jitter(width = 0.25), show.legend = FALSE, size = 3, alpha = 0.4) +
  geom_text(data = table_merg_f2, aes(x = as.factor(`input`), label = round((..count..)/2, digits = 0), y = (..count..)/2), stat= "count", hjust = 1.5, size = 5, color = "black", angle = 90) +
  ylab("quantified protein groups") + xlab("bead volume [無]") + theme_bw() +  scale_fill_manual(values = greens) +
  scale_fill_manual(values = greens) +
  theme(plot.title = element_text(color="black", size=16, face= "bold"), axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14), axis.title.y = element_text(color = "black", size = 14, face = "bold"))
quant_prot

accur_prot <- ggplot() +
  geom_boxplot(data = table_merg_f2, aes(x = as.factor(`input`), y = log2(`Ratio_norm_H/M`), fill = as.factor(`input`)), show.legend = FALSE, outlier.alpha = 0.15) +
  ylab("protein group\nlog2(SILAC H/M)") + xlab("bead volume [無]") + theme_bw() +  scale_fill_manual(values = greens) + ylim(-6,6) +
  theme(plot.title = element_text(color="black", size=16, face= "bold"), axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14), axis.title.y = element_text(color = "black", size = 14, face = "bold"))
accur_prot


### quant precision (CV)
prot_lf2 <- table_merg_f2

sd_l <- aggregate(data = prot_lf2, `Ratio_H/M` ~ `Gene_names` + `Protein_IDs` + `input`, FUN = sd, na.rm = FALSE) # select corret ratio!
colnames(sd_l) <- c("Gene_names","Protein.Group","sample","sd")

mean_l <- aggregate(data = prot_lf2, `Ratio_H/M` ~ `Gene_names` + `Protein_IDs` + `input`, FUN = mean, na.rm = FALSE) # select corret ratio!
colnames(mean_l) <- c("Gene_names","Protein.Group","sample","mean")

stats_l <- merge(mean_l, sd_l)
stats_l$cv <- (stats_l$sd / stats_l$mean)*100

p_meds <- aggregate(stats_l$`cv`, by = list(stats_l$`sample`), FUN = median, na.rm = TRUE)

CV <- ggplot(data = stats_l, aes(x = as.factor(`sample`), y = `cv`)) +
  geom_boxplot(aes(fill = as.factor(`sample`)), color = "black", show.legend = FALSE, outlier.alpha = 0.15) + 
  theme_bw() + xlab("bead volume [無]") + ylab("protein group CV [%]") + scale_y_continuous(limits = c(-12.5,150)) +
  scale_fill_manual(values = greens) +
  geom_text(data = p_meds, aes(x = as.factor(`Group.1`), y = -7.5, label = round(`x`, digits = 2)), size = 5, angle = 0) +
  theme(plot.title = element_text(color="black", size=12, face= "bold"), axis.title.x = element_text(color="black", size=14, face= "bold"),
        axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14, face="bold"))
CV

# create figure panel
quant_metrics_MAA <- plot_grid(SILAC_violin, quant_prot, accur_prot, CV, ncol = 2, align = "v")
ggdraw(quant_metrics_MAA)



#### analysis of data with protein input dilution series for automated enrichment
# subset evidence precursor level data to exclude MAA bead dilution data in this section
evidence_f <- evidence[grepl("ug_", evidence$Experiment),] # remove MAA bead dilution and non-enriched samples 
colnames(evidence_f)

# protein input dilution
neworder <- c("1","5","10","25","50","100","150","200","250","300")
evidence_f <- evidence_f %>% arrange(transform(evidence_f, `input` = factor(`input`, levels = neworder)), `input`)
evidence_f$input <- as.numeric(evidence_f$input)

ggplot(evidence_f, aes(x = `Retention_time`, fill = `Experiment`)) +
  geom_histogram(show.legend = FALSE, linewidth = 1, alpha = 0.2, bins = 105, color = "black") + #geom_smooth(show.legend = FALSE, color = "black", size = 1.5, alpha = 0.7) +
  scale_y_continuous(limits = c()) + facet_wrap('Experiment', ncol = 2) +
  theme_bw() + xlab("Retention time [min]") + ylab("identified precursors per min") +
  theme(axis.title.x = element_text(size = 14, face = "bold"), axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14, face = "bold"), axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold"))

# calculate S/N ratio
evidence_f$SILAC_enrich <- (evidence_f$`Ratio_H/L` + evidence_f$`Ratio_M/L`)
evidence_f$n <- rep(1)


# plot SILAC enrich ratios
p_meds <- aggregate(evidence_f$`SILAC_enrich`, by = list(evidence_f$`input`), FUN = median, na.rm = TRUE)

SILAC_violin <- ggplot(data = evidence_f, aes(x = as.factor(`input`), y = log2(`SILAC_enrich`))) +
  geom_hline(yintercept = 0, linetype="dashed", color = "black", size = 0.25) +
  geom_violin(aes(fill = as.factor(`input`)), show.legend = FALSE) + geom_boxplot(aes(fill = as.factor(`input`)), width = 0.2, show.legend = FALSE, outlier.size = -1) + 
  theme_bw() + xlab("protein input [痢]") + ylab("log2(labeled-/unlabeled\nprecursor intensity)") + scale_y_continuous(limits = c(min(log2(evidence_f$SILAC_enrich)), max(log2(evidence_f$SILAC_enrich), na.rm = T)+2)) +
  scale_fill_manual(values = blues) +
  geom_text(data = p_meds, aes(x = as.factor(`Group.1`), y = max(log2(evidence_f$SILAC_enrich), na.rm = T)+1, label = round(log2(`x`), digits = 2)), size = 5, vjust = 0.5, angle = 0) +
  theme(plot.title = element_text(color="black", size=12, face= "bold"), axis.title.x = element_text(color="black", size=14, face= "bold"),
        axis.text.x = element_text(color = "black", size=14, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14, face="bold"))
SILAC_violin

prec_HM <- evidence_f[,c("Experiment","input","n","Ratio_H/M")]
prec_HM[(prec_HM) == 0] <- NA
prec_HM <- prec_HM[!is.na(prec_HM$`Ratio_H/M`),] # remove NA values

all_precursors_agg <- aggregate(data = prec_HM, `n` ~ `Experiment` + `input`, FUN = sum, na.rm = TRUE)
all_precursors_agg$class <- rep("quantified")
all_precursors_agg$class <- rep("a")

p_meds <- aggregate(list(all_precursors_agg$`n`), by = list(all_precursors_agg$`input`), median)
colnames(p_meds) <- c("input","n")
p_meds$class <- rep("quantified")

quant_prec <- ggplot() +
  geom_path(data = all_precursors_agg, aes(x = `input`, y = `n`, color = as.factor(`class`), group = `class`), show.legend = FALSE, size = 1.5, alpha = 0.7) +
  geom_jitter(data = all_precursors_agg, aes(x = `input`, y = `n`, color = as.factor(`class`), shape = as.factor(`class`), group = `class`), position = position_jitter(width = 0.1), show.legend = TRUE, size = 3, alpha = 0.6) +
  geom_text_repel(data = p_meds, aes(x = `input`, label = round(n, digits = 0), y = `n`), vjust = 0.5, size = 5, color = "black", angle = 90) +
  ylab("number of precursors") + xlab("protein input [痢]") + theme_bw() + scale_shape_manual(values = c(19,17)) +
  scale_color_distiller(palette = "Blues", direction = 1) +theme(plot.title = element_text(color="black", size=16, face= "bold"), axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(color = "black", size = 14), legend.position = c(0.5,0.15), legend.background = element_blank(), legend.text = element_text(size = 14), legend.title = element_blank(),
        #axis.text.x = element_text(vjust = 0.5, color = "black", size=14, face = "bold", angle = 0, hjust = 0.5),
        axis.text.y = element_text(color = "black", size = 14), axis.title.y = element_text(color = "black", size = 14, face = "bold"))
quant_prec

ggplot() +
  geom_path(data = p_meds, aes(x = as.numeric(`input`), y = `n`, group = `class`), show.legend = FALSE, size = 1.5, alpha = 0.4, color = "black") +
  geom_point(data = all_precursors_agg, aes(x = as.numeric(`input`), y = `n`, color = `input`), show.legend = FALSE, size = 3, alpha = 1.0) +
  geom_point(data = all_precursors_agg, aes(x = as.numeric(`input`), y = `n`), show.legend = FALSE, size = 3, alpha = 0.8, shape = 1, color = "black") +
  geom_text_repel(data = p_meds, aes(x = as.numeric(`input`), label = round(`n`, digits = 0), y = round(`n`, digits = 0)), hjust = 1.5, size = 5, color = "black", angle = 90) +
  ylab("quantified precursors") + xlab("protein input [痢]") + theme_bw() +
  scale_color_distiller(palette = "Blues", direction = 1) + #ylim(c(0,3700)) +
  theme(plot.title = element_text(color="black", size=16, face= "bold"), axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        legend.position = c(0.5,0.15), legend.background = element_blank(), legend.text = element_text(size = 14), legend.title = element_blank(),
        axis.text.x = element_text(vjust = 0.5, color = "black", size=14, angle = 0, hjust = 0.5),
        axis.text.y = element_text(color = "black", size = 14), axis.title.y = element_text(color = "black", size = 14, face = "bold"))

## subset protein group, protein level data to exclude MAA bead dilution data in this section
table_merg_f <- read_tsv(file = "table_merg_f.tsv")
table_merg_f <- table_merg_f[grepl("ug_", table_merg_f$Experiment),] # filter for input dilution samples

# number of quantified proteins per sample
table_merg_f2 <- table_merg_f[!is.na(table_merg_f$`Ratio_norm_H/M`),] # remove 0 values for norm SILAC ratios

# 
prot_n <- aggregate(table_merg_f2$`Ratio_norm_H/M`, by = list(table_merg_f2$`Experiment`,table_merg_f2$`input`), FUN = length)
colnames(prot_n) <- c("Experiment","input","quantified_proteins")

# plot number of quantified proteins in bar chart
table_merg_f2$n <- rep(1)
all_prot_agg <- aggregate(data = table_merg_f2, `n` ~ `Experiment` + `input`, FUN = sum, na.rm = TRUE)

p_meds <- aggregate(list(all_prot_agg$`n`), by = list(all_prot_agg$`input`), median)
colnames(p_meds) <- c("input","n")
p_meds$group <- rep("quant")

quant_prot <- ggplot() +
  geom_path(data = p_meds, aes(x = as.numeric(`input`), y = `n`, group = `group`), show.legend = FALSE, size = 1.5, alpha = 0.4, color = "black") +
  geom_point(data = prot_n, aes(x = as.numeric(`input`), y = `quantified_proteins`, color = `input`), show.legend = FALSE, size = 3, alpha = 1.0) +
  geom_point(data = prot_n, aes(x = as.numeric(`input`), y = `quantified_proteins`), show.legend = FALSE, size = 3, alpha = 0.8, shape = 1, color = "black") +
  geom_text_repel(data = table_merg_f2, aes(x = as.numeric(`input`), label = round((..count..)/2, digits = 0), y = (..count..)/2), stat= "count", hjust = 1.5, size = 5, color = "black", angle = 90) +
  ylab("quantified protein groups") + xlab("protein input [痢]") + theme_bw() + #scale_shape_manual(values = c(19,17)) +
  scale_color_distiller(palette = "Blues", direction = 1) + ylim(c(0,3700)) +
  theme(plot.title = element_text(color="black", size=16, face= "bold"), axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        legend.position = c(0.5,0.15), legend.background = element_blank(), legend.text = element_text(size = 14), legend.title = element_blank(),
        axis.text.x = element_text(vjust = 0.5, color = "black", size=14, angle = 0, hjust = 0.5),
        axis.text.y = element_text(color = "black", size = 14), axis.title.y = element_text(color = "black", size = 14, face = "bold"))
quant_prot

accur_prot <- ggplot() +
  geom_boxplot(data = table_merg_f2, aes(x = as.factor(`input`), y = log2(`Ratio_norm_H/M`), fill = as.factor(`input`)), show.legend = FALSE, outlier.alpha = 0.15) +
  ylab("protein group\nlog2(SILAC H/M)") + xlab("protein input [痢]") + theme_bw() +  scale_fill_manual(values = blues) + ylim(-4,4) +
  theme(plot.title = element_text(color="black", size=16, face= "bold"), axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14), axis.title.y = element_text(color = "black", size = 14, face = "bold"))
accur_prot


### quant precision (CV)
prot_lf2 <- table_merg_f2

sd_l <- aggregate(data = prot_lf2, `Ratio_H/M` ~ `Gene_names` + `Protein_IDs` + `input`, FUN = sd, na.rm = FALSE) # select corret ratio!
colnames(sd_l) <- c("Gene_names","Protein.Group","sample","sd")

mean_l <- aggregate(data = prot_lf2, `Ratio_H/M` ~ `Gene_names` + `Protein_IDs` + `input`, FUN = mean, na.rm = FALSE) # select corret ratio!
colnames(mean_l) <- c("Gene_names","Protein.Group","sample","mean")

stats_l <- merge(mean_l, sd_l)
stats_l$cv <- (stats_l$sd / stats_l$mean)*100

p_meds <- aggregate(stats_l$`cv`, by = list(stats_l$`sample`), FUN = median, na.rm = TRUE)

CV <- ggplot(data = stats_l, aes(x = as.factor(`sample`), y = `cv`)) +
  geom_boxplot(aes(fill = as.factor(`sample`)), color = "black", show.legend = FALSE, outlier.alpha = 0.15) + 
  theme_bw() + xlab("protein input [痢]") + ylab("protein group CV [%]") + scale_y_continuous(limits = c(-25,100)) +
  scale_fill_manual(values = blues) +
  geom_text(data = p_meds, aes(x = as.factor(`Group.1`), y = -7.5, label = round(`x`, digits = 1)), size = 5, angle = 0) +
  theme(plot.title = element_text(color="black", size=12, face= "bold"), axis.title.x = element_text(color="black", size=14, face= "bold"),
        axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14, face="bold"))
CV

# create figure panel
quant_metrics_MAA <- plot_grid(SILAC_violin, quant_prot, accur_prot, CV, ncol = 2, align = "v")
ggdraw(quant_metrics_MAA)


