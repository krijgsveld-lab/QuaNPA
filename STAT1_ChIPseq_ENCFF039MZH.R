##### processing public STAT1 ChIPseq data from Hela cells treated with IFNg for 30 min #####
# data retrieved from ENCODE database - https://www.encodeproject.org/experiments/ENCSR000EZK/

### load necessary packages
library(BiocManager)
library(tidyverse)
library(Rsamtools)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(ggraph)
library(clusterProfiler)
library(ChIPpeakAnno)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(reshape2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(Biostrings)
library(ChIPseeker)
library(biomaRt)
library(phylotools)
library(AnnotationHub)

### set wd
setwd("~path to your working directory") # enter path here !

# load txdb gene model and shorten name for convenience
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

### load files (narrowpeak = BED format)
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")

# load bed file with conservative IDR thresholded peaks - https://www.encodeproject.org/files/ENCFF039MZH/
ChIP_bed <- import("ENCFF039MZH.bed", format = "BED", extraCols = extraCols_narrowPeak)

### use reference data for transcription start sites of human hg38 genome for annotation of peaks
data(TSS.human.GRCh38)
ucsc.hg38.knownGene <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

# annotate GRanges
ChIP.anno <- annotatePeakInBatch(ChIP_bed, AnnotationData = ucsc.hg38.knownGene)
ChIP.anno <- addGeneIDs(annotatedPeak = ChIP.anno, orgAnn = "org.Hs.eg.db", feature_id_type = "entrez_id", IDs2Add = c("ensembl","symbol"), mart = ensembl)

## calculate percentage of peaks in promoters etc.
ChIP_aCR <- assignChromosomeRegion(ChIP.anno, nucleotideLevel=FALSE, precedence=c("Promoters", "immediateDownstream", 
                                   "fiveUTRs", "threeUTRs", "Exons", "Introns"), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)


# create dataframes with percentages of peaks in different classes of genomic regios
ChIP_aCR_p <- ChIP_aCR$percentage
ChIP_aCR_p <- as.data.frame(ChIP_aCR_p)

# plot bargraphs
ggplot(ChIP_aCR_p, aes(x=`subjectHits`, y=`Freq`)) + geom_bar(position= "dodge", colour="black", width=.8, stat = "identity",
  show.legend = TRUE) + ylab("percentage [%]") + xlab("genomic elements") + theme_bw() + labs(title = "ChIPseq peaks") +
  theme(plot.title = element_text(color="black", size=16, face= "bold"),
        axis.title.x = element_text(color="black", size=16, face= "bold"),
        axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0.0, color="black", size=12, face= "bold"),
        axis.text.y = element_text(color="black", size=12, face = "bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"))

# peak heatmap
peakHeatmap <- peakHeatmap(ChIP.anno, TxDb = txdb, upstream = 2000, downstream = 2000, color = "firebrick3", title = "peak heatmap")

# density plot
avgprofile <- plotAvgProf2(ChIP.anno, TxDb = txdb, upstream = 2000, downstream = 2000, conf = 0.95, resample = 1000,
             xlab = "Genomic Region (5'->3')", ylab = "Read Count Frequency ChIPseq peaks")
avgprofile

### get promoters of hg38 ref genome via txdb object
genes <- genes(txdb)
promoters <- promoters(genes, upstream = 2000, downstream = 2000) #  define promoter range around TSS
promoters

# subset annotaded ChIPseq peaks within promoters
ChIP_prom <- subsetByOverlaps(ChIP.anno, promoters)
ChIP_prom

### create dataframe of GRanges object of ChIPseq data with peak sequence
ChIP_prom <- getAllPeakSequence(ChIP_prom, genome = Hsapiens)

ChIP_prom <- subsetByOverlaps(ChIP.anno, genes)
ChIP_prom <- getAllPeakSequence(ChIP.anno, genome = Hsapiens) # for unsubsetted data (not restricted to core promoter)

ChIP_prom_df <- as.data.frame(ChIP_prom)
ChIP_prom_df <- tibble::rownames_to_column(ChIP_prom_df, "Ranges")

# create fasta of peak seqeunces
ChIP_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, ChIP_prom)

# export fasta
writeXStringSet(ChIP_seq, "./ChIP_seq.fasta")

# import fasta sequences
ChIP_seq <- readDNAStringSet("ChIP_seq.fasta")

# create GC % values for each peak sequence
#count G and C nucleotides in sequence and divide by all nucleotides etc.
ChIP_prom_df_G <- str_count(ChIP_prom_df$sequence, "G")
ChIP_prom_df_C <- str_count(ChIP_prom_df$sequence, "C")
ChIP_prom_df_GCAT <- str_count(ChIP_prom_df$sequence, "")
ChIP_prom_df_GC <- ((ChIP_prom_df_C + ChIP_prom_df_G) / ChIP_prom_df_GCAT) * 100                                

# add % GC-content as column to dataframe
ChIP_prom_df["GC_content"] <- ChIP_prom_df_GC
colnames(ChIP_prom_df)

# save as tsv and import 
write_tsv(ChIP_prom_df, "ChIPseq_prom_df.tsv")
ChIP_prom_df <- read_tsv(file = "ChIPseq_prom_df.tsv")

# create genelist with 
STAT1_ChIPseq_targets <- ChIP_prom_df[,c("symbol","ensembl","signalValue")]

# save and import gene list
write_tsv(STAT1_ChIPseq_targets, "STAT1_ChIPseq_targets.tsv")
STAT1_ChIPseq_targets <- read_tsv(file = "STAT1_ChIPseq_targets.tsv")

no_dup_ChIP <- aggregate(data = STAT1_ChIPseq_targets, `signalValue` ~ `symbol`, FUN = max, na.rm = TRUE) # max signalValue for gene
no_dup_ChIP <- no_dup_ChIP[with(no_dup_ChIP, order(`signalValue`, decreasing = TRUE)),]
top_500_ChIP <- head(no_dup_ChIP, n = 500) # select top 500 STAT1 targets according to ChIP peak size

# save table of top 500 STAT1 ChIP targets
write_tsv(top_500_ChIP, "Hela_30min_IFNg_top_500_ChIP.tsv")
