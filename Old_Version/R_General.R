###This whole set of values state the whole variety of packages
###Will be used in order to become a complete analysis over the whole set.

###Inside of it we will see that the errors may ocurr but will be solved further on
library(readxl)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(tidyverse)
library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library("GenomicFeatures")
library(Rsubread)
library(UpSetR)
library(clusterProfiler) # Not working in the actual moment
library(biomaRt)
library(edgeR)
library(data.table)

getwd()
setwd("/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla")



RPKM_calc <- function(CountsXRead, Total_Mapped_Reads, Gene_len) {
  #Notes, it's important to have Gene len in bp
  num <- CountsXRead*10e9
  den <- Total_Mapped_Reads*Gene_len
  
}



###'
###'Now we will be importing the metadata, and the propper count files in order to do the analysis
###'

metadata <- read_excel("metadata_HEV_PHH.xlsx")

metadata <- metadata %>%
  dplyr::select("Run", "inoculation", "treatment", "time_point_of_rna_extraction_post_inoculation")

setwd(dir = "HTseqRes/HiSat2/")

New_counts <-list.files()[1:12]

for (file in New_counts) {
  data <- fread(file, header = FALSE, sep = "\t")
  
  parts <- strsplit(basename(file), "_")
  extracted_name <- parts[[1]][3] 
  

  assign(extracted_name, data)
  
}

SRR9937547<-SRR9937547 %>%
  rename( "SRR9937547" ="V2" )

SRR9937548<-SRR9937548 %>%
  rename( "SRR9937548" ="V2" )

SRR9937549<-SRR9937549 %>%
  rename( "SRR9937549" ="V2" )

SRR9937550<-SRR9937550 %>%
  rename( "SRR9937550" ="V2" )

SRR9937551<-SRR9937551 %>%
  rename( "SRR9937551" ="V2" )

SRR9937552<-SRR9937552 %>%
  rename( "SRR9937552" ="V2" )

SRR9937559<-SRR9937559 %>%
  rename( "SRR9937559" ="V2" )

SRR9937560<-SRR9937560 %>%
  rename( "SRR9937560" ="V2" )

SRR9937561<-SRR9937561 %>%
  rename( "SRR9937561" ="V2" )

SRR9937562<-SRR9937562 %>%
  rename( "SRR9937562" ="V2" )

SRR9937563<-SRR9937563 %>%
  rename( "SRR9937563" ="V2" )

SRR9937564<-SRR9937564 %>%
  rename( "SRR9937564" ="V2" )

###' Prepare the data table to work with, inside of it we will have the whole amount of information so we will
###' be able to perform the analysis further on in a much more automatically way, adapt the overall data and 
###' generate the metadata to work with DESeq2.

SAMPLES <- SRR9937547 %>% 
  dplyr::select(SRR9937547)%>%
  mutate(SRR9937548 = SRR9937548$SRR9937548, SRR9937549 = SRR9937549$SRR9937549,
         SRR9937550 = SRR9937550$SRR9937550, SRR9937551 = SRR9937551$SRR9937551,
         SRR9937552 = SRR9937552$SRR9937552, SRR9937559 = SRR9937559$SRR9937559,
         SRR9937560 = SRR9937560$SRR9937560, SRR9937561 = SRR9937561$SRR9937561,
         SRR9937562 = SRR9937562$SRR9937562, SRR9937563 = SRR9937563$SRR9937563,
         SRR9937564 = SRR9937564$SRR9937564)


SAMPLES <- data.frame(SAMPLES, row.names = SRR9937547$V1)
SAMPLES <- SAMPLES[1:(dim(SAMPLES)[1]-5),]


Used_metadata <- metadata %>% 
        filter(Run %in% colnames(SAMPLES))%>%
        mutate(inoculation =case_when(inoculation == "HEV genotype 3 p6 (Kernow C-1) G1634R MOI = 1" ~ "HEV_genotype_3_p6",
                                      inoculation == "control" ~ "control"))%>%
        mutate(across(where(is.character),as.factor)) 
  

###'Use the data prepared with DESeq and then perform the analysis of the data, consider
###'that the experimental design depending of the aim of the test. In this case we will look
###'after the differences between control and test
dsds <-DESeqDataSetFromMatrix(countData = SAMPLES,
                              colData = Used_metadata,
                              design = ~ inoculation )
dds <- DESeq(dsds)

results <- results(dds)
results_df<-as.data.frame(results)
results_df<-na.omit(results_df)

#filtered_results<- results_df[results_df$pvalue<0.05, ]
fc_filtered_res_down <- results_df[results_df$log2FoldChange <= -1, ]
fc_filtered_res_up <- results_df[results_df$log2FoldChange >= 1, ]


write.csv(results, file = "differential_expression_results.csv")


RPKM_calc <- function(CountsXRead, Total_Mapped_Reads, Gene_len) {
  num = CountsXRead*10e9
  den = Total_Mapped_Reads*Gene_len
  
  return(num/den)
  
}


GTF <- import.gff("../../HiSat2Res/Workfiles/genomic.gff", format="gff", feature = "exon")
GTF_exonic <- GTF[GTF$type == "exon", ]

exon_lengths <- data.frame(Gene = GTF_exonic@elementMetadata@listData[["gene"]], Length =  GTF@ranges@width)


exon_lengths<- exon_lengths %>%
  group_by(Gene) %>%
  summarise(Sum_Length = sum(Length))

exon_lengths <-data.frame(Length = exon_lengths$Sum_Length, row.names = exon_lengths$Gene)


##################

exon_lengths$Gene = rownames(exon_lengths)

merged_data <- merge(SAMPLES, exon_lengths, by.x = 0, by.y = "Gene", all.x = TRUE, all.y = TRUE)

rownames(merged_data) <- merged_data$Row.names

merged_data <- merged_data[, -1]

merged_data$Length[is.na(merged_data$Length)] <- 0

merged_data <- merged_data[rowSums(merged_data[])>0,]

merged_data <- na.omit(merged_data)

###'We may remove value 0 length? so we end up with capable values of meaningful length

#################
###'Working with rpkm values, calculating the mean of overall, individual and the control/infected group data.
###'Further on will allow us to play with the overall class of data.

sum_reads_per_gene <- rowSums(merged_data[, -ncol(merged_data)])

all_mapped_reads <- sum(sum_reads_per_gene)

rpkm_df <- data.frame(matrix(0, nrow = nrow(merged_data), ncol = ncol(merged_data)-1))
rownames(rpkm_df) <- rownames(merged_data)
colnames(rpkm_df) <- colnames(merged_data)[-ncol(merged_data)]

for (column in colnames(merged_data)[-ncol(merged_data)]) {
  sample_reads <- merged_data[,column]
  print(sample_reads)
  print(column)
  rpkm_df[, column] <- (sample_reads * 1e9) / (all_mapped_reads * merged_data$Length)
}

rpkm_df$Overall_RPKM <- rowMeans(rpkm_df[, -ncol(rpkm_df)])

rpkm_df <- rpkm_df[rowSums(rpkm_df[])>0,]

rpkm_df <- na.omit(rpkm_df)

rpkm_df <- rpkm_df[!is.infinite(rowSums(rpkm_df)),]

controlMeta <- Used_metadata[Used_metadata$inoculation == "control", ]

controlMeta$Run <- as.character(controlMeta$Run)

Control_RPKM<- rpkm_df[colnames(rpkm_df) %in% controlMeta$Run]
Infected_RPKM<- rpkm_df[!colnames(rpkm_df) %in% controlMeta$Run] 
Infected_RPKM <- Infected_RPKM[-ncol(Infected_RPKM)]

rpkm_df$Control_RPKM <- rowMeans(Control_RPKM)
rpkm_df$Infected_RPKM <- rowMeans(Infected_RPKM)

rpkm_means <- rpkm_df %>%
  dplyr::select(Overall_RPKM, Control_RPKM, Infected_RPKM)


rpkm_means <- rpkm_means[rowSums(rpkm_means[])>0,]

rpkm_means <- na.omit(rpkm_means)

rpkm_means <- rpkm_means[!is.infinite(rowSums(rpkm_means)),]

class(rpkm_means[1,1])


###'Filtering of the data, applying normalization and only keeping samples RPKM

sc <- log10(rpkm_df)+1
sc <- na.omit(sc)
sc <- sc[!is.infinite(rowSums(sc)),]
sc <- sc %>% filter(if_any(ends_with("_RPKM"),~ . >= 0.5)) 
sc <- sc[1:(length(sc)-3)]
sc1 <-  sc %>% top_n(1000)



###'Subset the first 4  and the 3 last timepoints in order to see the effect of the clustering information, so we can infere more 
###'how do they relate.
###'

Used_metadata <- Used_metadata %>%
  mutate(time_point_of_rna_extraction_post_inoculation = as.factor(str_replace_all(time_point_of_rna_extraction_post_inoculation, " ", "")))


First_Timepoints <- Used_metadata %>%
  filter(time_point_of_rna_extraction_post_inoculation %in% c("4h", "8h", "12h"))

Last_Timepoints <- Used_metadata %>%
  filter(!time_point_of_rna_extraction_post_inoculation %in% c("4h", "8h", "12h"))

rpkm_1st <- sc[First_Timepoints$Run]

rpkm_last <- sc[Last_Timepoints$Run]



###'To work with the different  fold changes into different plots we will select the different timepoints.
###'
###'
###

Timepoint1_SAMPLES<- SAMPLES[First_Timepoints$Run]

First_Timepoints$time_point_of_rna_extraction_post_inoculation <- droplevels(First_Timepoints$time_point_of_rna_extraction_post_inoculation)


dsds2_1st <-DESeqDataSetFromMatrix(countData = Timepoint1_SAMPLES,
                              colData = First_Timepoints,
                              design = ~ inoculation + time_point_of_rna_extraction_post_inoculation)
dds_1st <- DESeq(dsds2_1st)

results2 <- results(dds_1st)
results_df2<-as.data.frame(results2)
results_df2<-na.omit(results_df2)


Timepoint2_SAMPLES<- SAMPLES[Last_Timepoints$Run]


Last_Timepoints$time_point_of_rna_extraction_post_inoculation <- droplevels(Last_Timepoints$time_point_of_rna_extraction_post_inoculation)


dsds2_last <-DESeqDataSetFromMatrix(countData = Timepoint2_SAMPLES,
                                   colData = Last_Timepoints,
                                   design = ~ inoculation + time_point_of_rna_extraction_post_inoculation
                                   )
dds_last <- DESeq(dsds2_last)

results3 <- results(dds_last)
results_df3<-as.data.frame(results3)
results_df3<-na.omit(results_df3)


dsds3 <-DESeqDataSetFromMatrix(countData = SAMPLES,
                                    colData = Used_metadata,
                                    design = ~ inoculation + time_point_of_rna_extraction_post_inoculation
)
dds3 <- DESeq(dsds3)

results4 <- results(dds3)
results_df4<-as.data.frame(results4)
results_df4<-na.omit(results_df4)


BarplotSet <- data.frame(First_tp_FoldChange = results_df2$log2FoldChange)
rownames(BarplotSet) <- rownames(results_df2)

BarplotSet<- BarplotSet %>%
  mutate(Last_tp_FoldChange = ifelse(!is.na(match(rownames(BarplotSet), rownames(results_df3))), results_df3$log2FoldChange[match(rownames(BarplotSet), rownames(results_df3))], 0),
         All_tp_FoldChange = ifelse(!is.na(match(rownames(BarplotSet), rownames(results_df4))), results_df4$log2FoldChange[match(rownames(BarplotSet), rownames(results_df4),)], 0),
         GeneName = rownames(BarplotSet))

Long_BP <- tidyr::gather(BarplotSet, key = "TimePointSet", value = "FoldChange", -GeneName)


# Generate the count data estimated with vst so we can get the best propper meaning from the counts:

TP_1_vsd <- vst(dds_1st, blind = FALSE,nsub = 1000)
TP_2_vsd <- vst(dds_last, blind = FALSE,nsub = 1000)
TP_all_vsd <- vst(dds3, blind = FALSE,nsub = 1000)

TP_1_vsd_test <- as.data.frame(assay(TP_1_vsd)) %>% top_n(1000)
TP_2_vsd_test <- as.data.frame(assay(TP_2_vsd)) %>% top_n(1000)
TP_all_vsd_test <- as.data.frame(assay(TP_all_vsd)) %>% top_n(1000)

TP_1_vsd_nat <- as.data.frame(assay(TP_1_vsd)) 
TP_2_vsd_nat <- as.data.frame(assay(TP_2_vsd)) 
TP_all_vsd_nat <- as.data.frame(assay(TP_all_vsd)) 


###'====================================================================================================
###'                                        Plotting Block                                           ###
###'====================================================================================================
library(GGally)
library("EnhancedVolcano")

###'First we should determine which of the main thresholds it's important over the rest



pallete <- list(Inoculation = c( control = "#377bff", HEV_genotype_3_p6 = "#9e0034"), Time_Points = c("4h" = "#ff61b6", "8h" = "#53b847", "12h" = "#82138a", "24h" = "#787300", "48h" = "#af6ca4", "168h" = "#ff9066"))
sample_color_1 <- data.frame(Inoculation = Used_metadata$inoculation, Time_Points = Used_metadata$time_point_of_rna_extraction_post_inoculation)
rownames(sample_color_1) <- Used_metadata$Run 

pheatmap(as.matrix(sc), annotation_col =sample_color_1,cutree_cols = 2,main = "RPKM, top 100 genes untreated samples", annotation_colors = pallete ) 

##########

sample_color_1st <- data.frame(Inoculation = First_Timepoints$inoculation, Time_Points = First_Timepoints$time_point_of_rna_extraction_post_inoculation)
rownames(sample_color_1st) <- First_Timepoints$Run 

pheatmap(as.matrix(rpkm_1st), annotation_col =sample_color_1st,main = "RPKM,genes untreated samples 1st timepoints", annotation_colors = pallete ) 

sample_color_last <- data.frame(Inoculation = Last_Timepoints$inoculation, Time_Points = Last_Timepoints$time_point_of_rna_extraction_post_inoculation)
rownames(sample_color_last) <- Last_Timepoints$Run 

pheatmap(as.matrix(rpkm_last), annotation_col =sample_color_last,cutree_cols = 2,main = "RPKM,genes untreated samples last timepoints", annotation_colors = pallete ) 

###'LogFoldChange Barplot over timepoints.
ggplot(Long_BP, aes(x = GeneName, y = FoldChange, fill = TimePointSet)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("red", "lightgreen", "cadetblue1")) +
  labs(title = "Log2FoldChange, different timepoints",
       x = "Gene Name", y = "Fold Change") +
  ylim(-1.5,5)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggpairs(BarplotSet[-ncol(BarplotSet)])
ggpairs(results_df4)


pheatmap(as.matrix(sc), annotation_col =sample_color_1,cutree_cols = 2,main = "RPKM, top 100 genes untreated samples", annotation_colors = pallete ) 


pheatmap(as.matrix(TP_1_vsd_test), main = "First Timepoints Heatmap, normlized counts 1000 top samples",scale = "row", annotation_col = sample_color_1st,
         annotation_colors = pallete)
pheatmap(as.matrix(TP_2_vsd_test), main = "Last Timepoints Heatmap, normalized counts 1000 top samples", scale = "row", annotation_col = sample_color_last,
         annotation_colors = pallete)
pheatmap(as.matrix(TP_all_vsd_test), main = "All timepoints toghether Heatmap, normalized count 1000 top sampless", scale = "row", annotation_col = sample_color_1,
         annotation_colors = pallete, cluster_rows = FALSE)

pheatmap(as.matrix(TP_1_vsd_nat), main = "First Timepoints Heatmap, normlized counts s",scale = "row", annotation_col = sample_color_1st,
         annotation_colors = pallete)
pheatmap(TP_2_vsd_nat, main = "Last Timepoints Heatmap, normalized counts ", scale = "row", annotation_col = sample_color_last,
         annotation_colors = pallete)
pheatmap(TP_all_vsd_nat, main = "All timepoints toghether Heatmap, normalized count ", scale = "row", annotation_col = sample_color_1,
         annotation_colors = pallete)

# pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))
# ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#   coord_fixed()

###'============================================================================================'###
###'                                  STAR SAMPLING                                             '###
###'============================================================================================'###


###'Introduce the data of STAR inside the R part.
###'Then perform the hole analysis again.
###'

setwd(dir = "/mnt/Viro_Data/Mitarbeiter/Ian/Run_2_PNAS_paper/HTseqRes/STAR")

New_counts <-list.files()[1:12]

for (file in New_counts) {
  data <- fread(file, header = FALSE, sep = "\t")
  
  parts <- strsplit(basename(file), "_")
  extracted_name <- parts[[1]][3] 
  
  
  assign(extracted_name, data)
  
}

SRR9937547<-SRR9937547 %>%
  rename( "SRR9937547" ="V2" )

SRR9937548<-SRR9937548 %>%
  rename( "SRR9937548" ="V2" )

SRR9937549<-SRR9937549 %>%
  rename( "SRR9937549" ="V2" )

SRR9937550<-SRR9937550 %>%
  rename( "SRR9937550" ="V2" )

SRR9937551<-SRR9937551 %>%
  rename( "SRR9937551" ="V2" )

SRR9937552<-SRR9937552 %>%
  rename( "SRR9937552" ="V2" )

SRR9937559<-SRR9937559 %>%
  rename( "SRR9937559" ="V2" )

SRR9937560<-SRR9937560 %>%
  rename( "SRR9937560" ="V2" )

SRR9937561<-SRR9937561 %>%
  rename( "SRR9937561" ="V2" )

SRR9937562<-SRR9937562 %>%
  rename( "SRR9937562" ="V2" )

SRR9937563<-SRR9937563 %>%
  rename( "SRR9937563" ="V2" )

SRR9937564<-SRR9937564 %>%
  rename( "SRR9937564" ="V2" )

###' Prepare the data table to work with, inside of it we will have the whole amount of information so we will
###' be able to perform the analysis further on in a much more automatically way, adapt the overall data and 
###' generate the metadata to work with DESeq2.

SAMPLES_STAR <- SRR9937547 %>% 
  dplyr::select(SRR9937547)%>%
  mutate(SRR9937548 = SRR9937548$SRR9937548, SRR9937549 = SRR9937549$SRR9937549,
         SRR9937550 = SRR9937550$SRR9937550, SRR9937551 = SRR9937551$SRR9937551,
         SRR9937552 = SRR9937552$SRR9937552, SRR9937559 = SRR9937559$SRR9937559,
         SRR9937560 = SRR9937560$SRR9937560, SRR9937561 = SRR9937561$SRR9937561,
         SRR9937562 = SRR9937562$SRR9937562, SRR9937563 = SRR9937563$SRR9937563,
         SRR9937564 = SRR9937564$SRR9937564)


SAMPLES_STAR <- data.frame(SAMPLES_STAR, row.names = SRR9937547$V1)
SAMPLES_STAR <- SAMPLES_STAR[1:(dim(SAMPLES_STAR)[1]-5),]


###'============================================================================================'###
###'                                  BOWTIE SAMPLING                                           '###
###'============================================================================================'###


setwd(dir = "../BOWTIE/")

New_counts <-list.files()[2:13]

for (file in New_counts) {
  data <- fread(file, header = FALSE, sep = "\t")
  
  parts <- strsplit(basename(file), "_")
  extracted_name <- parts[[1]][3] 
  
  
  assign(extracted_name, data)
  
}

SRR9937547<-SRR9937547 %>%
  rename( "SRR9937547" ="V2" )

SRR9937548<-SRR9937548 %>%
  rename( "SRR9937548" ="V2" )

SRR9937549<-SRR9937549 %>%
  rename( "SRR9937549" ="V2" )

SRR9937550<-SRR9937550 %>%
  rename( "SRR9937550" ="V2" )

SRR9937551<-SRR9937551 %>%
  rename( "SRR9937551" ="V2" )

SRR9937552<-SRR9937552 %>%
  rename( "SRR9937552" ="V2" )

SRR9937559<-SRR9937559 %>%
  rename( "SRR9937559" ="V2" )

SRR9937560<-SRR9937560 %>%
  rename( "SRR9937560" ="V2" )

SRR9937561<-SRR9937561 %>%
  rename( "SRR9937561" ="V2" )

SRR9937562<-SRR9937562 %>%
  rename( "SRR9937562" ="V2" )

SRR9937563<-SRR9937563 %>%
  rename( "SRR9937563" ="V2" )

SRR9937564<-SRR9937564 %>%
  rename( "SRR9937564" ="V2" )

###' Prepare the data table to work with, inside of it we will have the whole amount of information so we will
###' be able to perform the analysis further on in a much more automatically way, adapt the overall data and 
###' generate the metadata to work with DESeq2.

SAMPLES_BOWTIE <- SRR9937547 %>% 
  dplyr::select(SRR9937547)%>%
  mutate(SRR9937548 = SRR9937548$SRR9937548, SRR9937549 = SRR9937549$SRR9937549,
         SRR9937550 = SRR9937550$SRR9937550, SRR9937551 = SRR9937551$SRR9937551,
         SRR9937552 = SRR9937552$SRR9937552, SRR9937559 = SRR9937559$SRR9937559,
         SRR9937560 = SRR9937560$SRR9937560, SRR9937561 = SRR9937561$SRR9937561,
         SRR9937562 = SRR9937562$SRR9937562, SRR9937563 = SRR9937563$SRR9937563,
         SRR9937564 = SRR9937564$SRR9937564)


SAMPLES_BOWTIE <- data.frame(SAMPLES_BOWTIE, row.names = SRR9937547$V1)
SAMPLES_BOWTIE <- SAMPLES_BOWTIE[1:(dim(SAMPLES_BOWTIE)[1]-5),]

###' Plots for presentations:

SAMPLING_BW <- SAMPLES_BOWTIE
SAMPLING <- SAMPLES
SAMPLING_STAR <- SAMPLES_STAR



SAMPLING_BW <- SAMPLING_BW[!rowSums(SAMPLING_BW == 0) == ncol(SAMPLING_BW), ]
SAMPLING <- SAMPLING[!rowSums(SAMPLING == 0) == ncol(SAMPLING), ]
SAMPLING_STAR <- SAMPLING_STAR[!rowSums(SAMPLING_STAR == 0) == ncol(SAMPLING_STAR), ]


SAMPLING_BW$Gene <- rownames(SAMPLING_BW)
SAMPLING$Gene <- rownames(SAMPLING)
SAMPLING_STAR$Gene <- rownames(SAMPLING_STAR)


BW <- gather(SAMPLING_BW, key = "Sample", value = "Count", -Gene)
STAR <- gather(SAMPLING_STAR, key = "Sample", value = "Count", -Gene)
HSAT <- gather(SAMPLING, key = "Sample", value = "Count", -Gene)

BW <- BW %>% 
  dplyr::arrange(desc(Count)) %>%
  head(30000)

STAR <- STAR %>% 
  dplyr::arrange(desc(Count)) %>%
  head(30000)

HSAT <- HSAT %>% 
  dplyr::arrange(desc(Count)) %>%
  head(30000)


ggplot() + 
  geom_point(data = BW, aes(y = Count, x = reorder(Gene, Count), shape = Sample, color = "Bowtie2"))+
  geom_point(data = HSAT, aes(y = Count, x = reorder(Gene, Count), shape = Sample, color = "HiSat2" ))+
  geom_point(data = STAR, aes(y = Count, x = reorder(Gene, Count), shape = Sample, color = "STAR" ))+
  scale_shape_manual(values = 1:12) + scale_color_manual(values = c( "#ff6e95", "#bc49b7", "#a9c93d"))+
  labs(title = "Count comparison between different gene mapping tools, 1000 top",
       x = "Gene",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom",
        legend.title = element_blank(), 
        axis.text.x = element_blank())


ggplot() + 
  geom_line(data = BW, aes(x = Count, y = rank(Count), color = "Bowtie2"))+
  geom_line(data = HSAT, aes(x = Count, y = rank(Count), color = "HiSat2" ))+
  geom_line(data = STAR, aes(x = Count, y = rank(Count), color = "STAR" ))+ 
  geom_point(data = BW, aes(x = Count, y =  rank(Count), shape = Sample, color = "Bowtie2"))+
  geom_point(data = HSAT, aes(x = Count, y =  rank(Count), shape = Sample, color = "HiSat2" ))+
  geom_point(data = STAR, aes(x = Count, y =  rank(Count), shape = Sample, color = "STAR" ))+
  # geom_smooth(method = "lm", data = BW, aes(x = Count, y = rank(Count), color = "Bowtie2")) + 
  # geom_smooth(method = "lm", data = HSAT, aes(x = Count, y = rank(Count), color = "HiSat2" ))+
  # geom_smooth(method = "lm", data = STAR, aes(x = Count, y = rank(Count), color = "STAR" ))+
  scale_color_manual(values = c( "#ff6e95", "#bc49b7", "#a9c93d"))+ scale_shape_manual(values = 1:12)+ 
  labs(title = "Count comparison between different gene mapping tools",
       x = "Count",
       y = "Rank") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom",
        legend.title = element_blank())



CorrelationSet <- data.frame(HiSat_SRR9937550 = SAMPLING$SRR9937550, HiSat_SRR9937548 = SAMPLING$SRR9937548 , HiSat_SRR9937562 = SAMPLING$SRR9937562,  HiSat_SRR9937563 = SAMPLING$SRR9937563)
rownames(CorrelationSet) <- rownames(SAMPLING) 

CorrelationSet<- CorrelationSet %>% 
mutate(STAR_SRR9937550 = ifelse(!is.na(match(rownames(CorrelationSet), rownames(SAMPLES_STAR))), SAMPLES_STAR$SRR9937550[match(rownames(CorrelationSet), rownames(SAMPLES_STAR))], 0),
       STAR_SRR9937548 = ifelse(!is.na(match(rownames(CorrelationSet), rownames(SAMPLES_STAR))), SAMPLES_STAR$SRR9937548[match(rownames(CorrelationSet), rownames(SAMPLES_STAR))], 0),
       STAR_SRR9937562 = ifelse(!is.na(match(rownames(CorrelationSet), rownames(SAMPLES_STAR))), SAMPLES_STAR$SRR9937562[match(rownames(CorrelationSet), rownames(SAMPLES_STAR))], 0),
       STAR_SRR9937563 = ifelse(!is.na(match(rownames(CorrelationSet), rownames(SAMPLES_STAR))), SAMPLES_STAR$SRR9937563[match(rownames(CorrelationSet), rownames(SAMPLES_STAR))], 0),
       Bowtie_SRR9937550 = ifelse(!is.na(match(rownames(CorrelationSet), rownames(SAMPLES_BOWTIE))), SAMPLES_BOWTIE$SRR9937550[match(rownames(CorrelationSet), rownames(SAMPLES_BOWTIE))], 0),
       Bowtie_SRR9937548 = ifelse(!is.na(match(rownames(CorrelationSet), rownames(SAMPLES_BOWTIE))), SAMPLES_BOWTIE$SRR9937548[match(rownames(CorrelationSet), rownames(SAMPLES_BOWTIE))], 0),
       Bowtie_SRR9937562 = ifelse(!is.na(match(rownames(CorrelationSet), rownames(SAMPLES_BOWTIE))), SAMPLES_BOWTIE$SRR9937562[match(rownames(CorrelationSet), rownames(SAMPLES_BOWTIE))], 0),
       Bowtie_SRR9937563 = ifelse(!is.na(match(rownames(CorrelationSet), rownames(SAMPLES_BOWTIE))), SAMPLES_BOWTIE$SRR9937563[match(rownames(CorrelationSet), rownames(SAMPLES_BOWTIE))], 0))

CorrelationSet <- CorrelationSet[rowSums(CorrelationSet) != 0, ]
CorrelationSet$Best_Tool <- colnames(CorrelationSet)[apply(CorrelationSet, 1, which.max)]

ggpairs(CorrelationSet, columns = 1:12, mapping = aes(color = Best_Tool)) 

###' Use this command in order to obtain random rows(genes) of the desired DF
###' We will generate the same rowplots as previously done so we will get the best position
###' of X random genes so we can see and recognise how do they vary.


RandomSelection <- sample(rownames(SAMPLING), 10)

RandomSelectionDF_BW <- BW %>%
  filter(Gene %in% RandomSelection)
RandomSelectionDF_STAR <- STAR %>%
  filter(Gene %in% RandomSelection)
RandomSelectionDF_HSAT <- HSAT %>%
  filter(Gene %in% RandomSelection)

ggplot() + 
  geom_point(data = RandomSelectionDF_BW, aes(y = Count, x = reorder(Gene, Count), shape = Sample, color = "Bowtie2"))+
  geom_point(data = RandomSelectionDF_HSAT, aes(y = Count, x = reorder(Gene, Count), shape = Sample, size =  color = "HiSat2" ))+
  geom_point(data = RandomSelectionDF_STAR, aes(y = Count, x = reorder(Gene, Count), shape = Sample, color = "STAR" ))+
  scale_shape_manual(values = 1:12) + scale_color_manual(values = c( "#ff6e95", "#bc49b7", "#a9c93d"))+
  labs(title = "Count comparison between different gene mapping tools, subsamples",
       x = "Gene",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom",
        legend.title = element_blank(), 
        axis.text.x = element_text(angle = 60, vjust = 0.7, hjust = 1))

###' We select the top, the lowest and  4 at random
Pseudorandom <- c(rownames(top_n(unique(SAMPLING_BW), 2)),sample(rownames(SAMPLING_BW), 4), rownames(tail(SAMPLES_BOWTIE, 2)))
PseudoRandomSelectionDF_BW <- BW %>%
  filter(Gene %in% Pseudorandom)
PseudoRandomSelectionDF_STAR <- STAR %>%
  filter(Gene %in% Pseudorandom)
PseudoRandomSelectionDF_HSAT <- HSAT %>%
  filter(Gene %in% Pseudorandom)

ggplot() + 
  geom_point(data = PseudoRandomSelectionDF_BW, aes(y = Count, x = reorder(Gene, Count), shape = Sample, size = "Bowtie2", color = "Bowtie2"))+
  geom_point(data = PseudoRandomSelectionDF_HSAT, aes(y = Count, x = reorder(Gene, Count), shape = Sample, size = "HiSat2", color = "HiSat2" ))+
  geom_point(data = PseudoRandomSelectionDF_STAR, aes(y = Count, x = reorder(Gene, Count), shape = Sample, size = "STAR", color = "STAR" ))+
  scale_shape_manual(values = 1:12) + scale_color_manual(values = c( "#ff6e95", "#bc49b7", "#a9c93d"))+ 
  labs(title = "Count comparison between different gene mapping tools, Pseudo-Random subsamples",
       x = "Gene",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom",
        legend.title = element_blank(), 
        axis.text.x = element_text(angle = 60, vjust = 0.7, hjust = 1))

