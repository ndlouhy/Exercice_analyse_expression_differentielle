elibrary("DESeq2")
library("ggplot2")
library("NMF")
library("ggbeeswarm")
library("genefilter")
library("pheatmap")
library(ggrepel)
library(edgeR)
library(gplots)
library(RColorBrewer)

### Open file

file_path <- "C:/Users/noedl/Documents/cours_M_bioinformatique/cours_master_bioinformatique/M3/recherche_stage_2025_26/institut_thorax_igdr/analyse_diff/data/gene_count.xls"

data_file <- read.table(file_path, header = T, sep = "\t")

### File formatting
rownames(data_file) <- data_file$gene_id
data_file <- data_file[,-1]


# ReadCount creation
ReadCount <- data_file[, grep("^Sample_", colnames(data_file))]
dim(ReadCount)

# Metadata creation
meta <- data.frame(
  sample = colnames(ReadCount),
  group = c(
    rep("Group1", 8),
    rep("Group2", 9),
    rep("Group3", 9)
  )
)

# Gene description file creation
genes_description <- data_file[, !grepl("^Sample_", colnames(data_file))]

## Data filtering
#Comptage des gènes non exprimés dans aucun échantillon
table(rowSums(ReadCount) == 0)

# Suppression de ces gènes
ReadCount <- ReadCount[rowSums(ReadCount) > 0, ]

# Vérification
table(rowSums(ReadCount) == 0)

cutoff <- cpm(10, mean(colSums(ReadCount)))
dim(ReadCount)
keep <- rowSums(cpm(ReadCount)>cutoff[1]) >= 2
summary(keep)
ReadCount <- ReadCount[keep,]


## Pre-analysis steps
# DESeq object creation
DESeq.ds <- DESeqDataSetFromMatrix (countData = ReadCount,
                                    colData = meta,
                                    design = ~ group)

colSums(counts(DESeq.ds))

# Normalization
DESeq.ds <- estimateSizeFactors(DESeq.ds)
DESeq.ds@colData$sizeFactor

# Log transformation
rld <- rlog(DESeq.ds, blind = FALSE)
rlog.norm.counts <- assay(rld)


# Histogram
distance.m_rlog  <- as.dist(1 - cor(rlog.norm.counts , method = "pearson" ))
plot(hclust(distance.m_rlog), labels = colnames(rlog.norm.counts),
     main = "Dendrogramme représentant la distance entre les échantillons\ndistance: Pearson  correlation")

# PCA
pcaData <- plotPCA(rld, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size=3) +
  geom_text_repel(size=3, max.overlaps=10) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA des échantillons par regroupement par groupe expérimental")


# Samples 2, 3, 7 and 12 appear to be outliers
# Samples suppression :
outliers <- c("Sample_2", "Sample_3", "Sample_7", "Sample_12")

ReadCount_filtered <- ReadCount[, !(colnames(ReadCount) %in% outliers)]
meta_filtered <- meta[!(meta$sample %in% outliers), ]

# Verification
DESeq.ds <- DESeqDataSetFromMatrix (countData = ReadCount_filtered,
                                    colData = meta_filtered,
                                    design = ~ group)

colSums(counts(DESeq.ds))

DESeq.ds <- estimateSizeFactors(DESeq.ds)
DESeq.ds@colData$sizeFactor


rld <- rlog(DESeq.ds, blind = FALSE)
rlog.norm.counts <- assay(rld)

#pdf("C:/Users/noedl/Documents/cours_M_bioinformatique/cours_master_bioinformatique/M3/recherche_stage_2025_26/institut_thorax_igdr/analyse_diff/Plots_results.pdf")

# Histogram
distance.m_rlog  <- as.dist(1 - cor(rlog.norm.counts , method = "pearson" ))
plot(hclust(distance.m_rlog), labels = colnames(rlog.norm.counts),
     main = "Dendrogramme représentant la distance entre les échantillons (sans outliers)\ndistance: Pearson  correlation")

# PCA
pcaData <- plotPCA(rld, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size=3) +
  geom_text_repel(size=3, max.overlaps=10) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA des échantillons par regroupement par groupe expérimental (sans outliers)")

#-----------------------------------------------------------------------
#============   Differential Expression with DESeq2     ================
#-----------------------------------------------------------------------

reference_group <- "Group1"
DESeq.ds$group <- factor(DESeq.ds$group, levels = c(reference_group, setdiff(levels(DESeq.ds$group), reference_group)))

DESeq.ds <- DESeq(DESeq.ds)
plotDispEsts(DESeq.ds)


group1 <- "Group1"
group2 <- "Group2"
group3 <- "Group3"


# Group 1 vs Group 2
DGE.results_1v2 <- results(DESeq.ds, c("group", group1, group2), pAdjustMethod = "BH")
summary(DGE.results_1v2)
mcols(DGE.results_1v2, use.names = TRUE)

# Group 1 vs Group 3
DGE.results_1v3 <- results(DESeq.ds, c("group", group1, group3), pAdjustMethod = "BH")
summary(DGE.results_1v3)
mcols(DGE.results_1v3, use.names = TRUE)

# Group 2 vs Group 3
DGE.results_2v3 <- results(DESeq.ds, c("group", group2, group3), pAdjustMethod = "BH")
summary(DGE.results_2v3)
mcols(DGE.results_2v3, use.names = TRUE)

# results
results.DESeq2_1v2 <- DGE.results_1v2
results.DESeq2_1v3 <- DGE.results_1v3
results.DESeq2_2v3 <- DGE.results_2v3


# MA Plot
# 1v2
DESeq2::plotMA(results.DESeq2_1v2, ylim = c(-4,4), main="MA plot entre le groupe 1 et le  groupe 2", alpha=0.0001)
abline(h = c(-1,1), col = "blue", lty = 2)
mtext(c("-1 fold", "+1 fold"), side = 4, at = c(-1, 1), cex = 0.8, line = 0.5, col = "blue")
legend("topright", legend=c("Gène", "Gène diff. exprimé"), col=c("grey", "blue"), pch=20)


# 1v3
DESeq2::plotMA(results.DESeq2_1v3, ylim = c(-4,4), main="MA plot entre le groupe 1 et le  groupe 3", alpha=0.0001)
abline(h = c(-1,1), col = "blue", lty = 2)
mtext(c("-1 fold", "+1 fold"), side = 4, at = c(-1, 1), cex = 0.8, line = 0.5, col = "blue")
legend("topright", legend=c("Gène", "Gène diff. exprimé"), col=c("grey", "blue"), pch=20)


# 2v3
DESeq2::plotMA(results.DESeq2_2v3, ylim = c(-4,4), main="MA plot: Group2 vs Group3", alpha=0.0001)
abline(h = c(-1,1), col = "blue", lty = 2)
mtext(c("-1 fold", "+1 fold"), side = 4, at = c(-1, 1), cex = 0.8, line = 0.5, col = "blue")
legend("topright", legend=c("Gène", "Gène diff. exprimé"), col=c("grey", "blue"), pch=20)



# Read counts of the top gene and plot, mainly useful that there was no mistake when setting the constrasts
topGene1v2 <- rownames(DGE.results_1v2)[which.min(DGE.results_1v2$padj)]
plotCounts(DESeq.ds, gene = topGene1v2, intgroup="group")

topGene1v3 <- rownames(DGE.results_1v3)[which.min(DGE.results_1v3$padj)]
plotCounts(DESeq.ds, gene = topGene1v3, intgroup="group")

topGene2v3 <- rownames(DGE.results_2v3)[which.min(DGE.results_2v3$padj)]
plotCounts(DESeq.ds, gene = topGene2v3, intgroup="group")


### Gene clustering, heatmap

# R-log transformation of the raw read counts
rld <- rlog(DESeq.ds, blind = FALSE)

alpharisk <- 0.001

# 1v2
DGE.results.sorted_1v2  <- DGE.results_1v2[order(DGE.results_1v2$padj), ]
# identify  genes  with  the  desired  adjusted p-value cut -off
DGEgenes_1v2  <- rownames(subset(DGE.results.sorted_1v2, padj < alpharisk))

########### RESULTS ###########

DESeq_results_1v2 <- DGE.results_1v2[DGE.results_1v2$padj < alpharisk & !is.na(DGE.results_1v2$padj),] #liste des gènes correspondant aux DE significatifs

DESeq_results_sorted_1v2 <- DESeq_results_1v2[order(DESeq_results_1v2$padj),]

DESeq_table_1v2 <- data.frame(
  Geneid = rownames(DESeq_results_sorted_1v2),
  GeneName = genes_description[rownames(DESeq_results_sorted_1v2), "gene_name"],
  baseMean = DESeq_results_sorted_1v2$baseMean,
  log2FoldChange = DESeq_results_sorted_1v2$log2FoldChange,
  lfcSE = DESeq_results_sorted_1v2$lfcSE,
  stat = DESeq_results_sorted_1v2$stat,
  pvalue = DESeq_results_sorted_1v2$pvalue,
  padj = DESeq_results_sorted_1v2$padj
)

# Creation of a descriptive table of differentially expressed genes
DESeq_table_1v2_annotated <- cbind(
  DESeq_table_1v2,
  genes_description[match(DESeq_table_1v2$Geneid, rownames(genes_description)), ]
)

# Check the result
head(DESeq_table_1v2_annotated)

# 1v3
DGE.results.sorted_1v3  <- DGE.results_1v3[order(DGE.results_1v3$padj), ]
# identify  genes  with  the  desired  adjusted p-value cut -off
DGEgenes_1v3  <- rownames(subset(DGE.results.sorted_1v3, padj < alpharisk))

########### RESULTS ###########

DESeq_results_1v3 <- DGE.results_1v3[DGE.results_1v3$padj < alpharisk & !is.na(DGE.results_1v3$padj),] #liste des gènes correspondant aux DE significatifs

DESeq_results_sorted_1v3 <- DESeq_results_1v3[order(DESeq_results_1v3$padj),]

DESeq_table_1v3 <- data.frame(
  Geneid = rownames(DESeq_results_sorted_1v3),
  GeneName = genes_description[rownames(DESeq_results_sorted_1v3), "gene_name"],
  baseMean = DESeq_results_sorted_1v3$baseMean,
  log2FoldChange = DESeq_results_sorted_1v3$log2FoldChange,
  lfcSE = DESeq_results_sorted_1v3$lfcSE,
  stat = DESeq_results_sorted_1v3$stat,
  pvalue = DESeq_results_sorted_1v3$pvalue,
  padj = DESeq_results_sorted_1v3$padj
)

# Creation of a descriptive table of differentially expressed genes
DESeq_table_1v3_annotated <- cbind(
  DESeq_table_1v3,
  genes_description[match(DESeq_table_1v3$Geneid, rownames(genes_description)), ]
)

# Check the result
head(DESeq_table_1v3_annotated)

# 2v3
DGE.results.sorted_2v3  <- DGE.results_2v3[order(DGE.results_2v3$padj), ]
# identify  genes  with  the  desired  adjusted p-value cut -off
DGEgenes_2v3  <- rownames(subset(DGE.results.sorted_2v3, padj < alpharisk))

########### RESULTS ###########

DESeq_results_2v3 <- DGE.results_2v3[DGE.results_2v3$padj < alpharisk & !is.na(DGE.results_2v3$padj),] #liste des gènes correspondant aux DE significatifs

DESeq_results_sorted_2v3 <- DESeq_results_2v3[order(DESeq_results_2v3$padj),]

DESeq_table_2v3 <- data.frame(
  Geneid = rownames(DESeq_results_sorted_2v3),
  GeneName = genes_description[rownames(DESeq_results_sorted_2v3), "gene_name"],
  baseMean = DESeq_results_sorted_2v3$baseMean,
  log2FoldChange = DESeq_results_sorted_2v3$log2FoldChange,
  lfcSE = DESeq_results_sorted_2v3$lfcSE,
  stat = DESeq_results_sorted_2v3$stat,
  pvalue = DESeq_results_sorted_2v3$pvalue,
  padj = DESeq_results_sorted_2v3$padj
)

# Creation of a descriptive table of differentially expressed genes
DESeq_table_2v3_annotated <- cbind(
  DESeq_table_2v3,
  genes_description[match(DESeq_table_2v3$Geneid, rownames(genes_description)), ]
)

# Check the result
head(DESeq_table_2v3_annotated)

## Heatmaps

# Extracting rlog data
rlog_mat <- assay(rld)

# Color 
color_palette <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)

# Extract sample metadata from DESeq object
coldata <- as.data.frame(colData(DESeq.ds))
# Selection of samples from group 1 and 2
samples_group_1v2 <- coldata$group %in% c("Group1", "Group2")
rlog_mat_1v2 <- rlog_mat[, samples_group_1v2]

# Extraction of differentially expressed genes between Group 1 and 2
DE_rlog_data_1v2 <- rlog_mat_1v2[intersect(rownames(rlog_mat_1v2), DGEgenes_1v2), ]
rownames(DE_rlog_data_1v2) <- genes_description[rownames(DE_rlog_data_1v2), "gene_name"]
DE_rlog_matrix_1v2 <- as.matrix(DE_rlog_data_1v2)

# Heatmap generation
if (nrow(DE_rlog_matrix_1v2) > 1) {
  heatmap.2(
    DE_rlog_matrix_1v2, 
    scale = "row",
    trace = "none",
    col = color_palette,
    main = "Heatmap entre le groupe 1 et le groupe 2",
    
    # Personalization
    labRow = NA,
    cexCol = 1,        
    srtCol = 45,  
    margins = c(8, 10)
  )
} else {
  cat("Aucun gène significatif pour Groupe 1 vs Groupe 2\n")
}

#-------------------------------------------
#  HEATMAP pour comparaison Group1 vs Group3
#-------------------------------------------
# Extract sample metadata from DESeq object
coldata <- as.data.frame(colData(DESeq.ds))
# Selection of samples from group 1 and 3
samples_group_1v3 <- coldata$group %in% c("Group1", "Group3")
rlog_mat_1v3 <- rlog_mat[, samples_group_1v3]

# Extraction of differentially expressed genes between Group 1 and 3
DE_rlog_data_1v3 <- rlog_mat_1v3[intersect(rownames(rlog_mat_1v3), DGEgenes_1v3), ]
rownames(DE_rlog_data_1v3) <- genes_description[rownames(DE_rlog_data_1v3), "gene_name"]
DE_rlog_matrix_1v3 <- as.matrix(DE_rlog_data_1v3)

# Heatmap generation
if (nrow(DE_rlog_matrix_1v3) > 1) {
  heatmap.2(
    DE_rlog_matrix_1v3, 
    scale = "row",
    trace = "none",
    col = color_palette,
    main = "Heatmap entre le groupe 1 et le groupe 3",
    
    # Personalization
    labRow = NA,
    cexCol = 1,
    srtCol = 45,  
    margins = c(8, 10) 
  )
} else {
  cat("Aucun gène significatif pour Groupe 1 vs Groupe 3\n")
}

#-------------------------------------------
#  HEATMAP pour comparaison Group2 vs Group3
#-------------------------------------------
# Extract sample metadata from DESeq object
coldata <- as.data.frame(colData(DESeq.ds))
# Selection of samples from group 2 and 3
samples_group_2v3 <- coldata$group %in% c("Group2", "Group3")
rlog_mat_2v3 <- rlog_mat[, samples_group_2v3]

# Extraction of differentially expressed genes between Group 2 and 3
DE_rlog_data_2v3 <- rlog_mat_2v3[intersect(rownames(rlog_mat_2v3), DGEgenes_2v3), ]
rownames(DE_rlog_data_2v3) <- genes_description[rownames(DE_rlog_data_2v3), "gene_name"]
DE_rlog_matrix_2v3 <- as.matrix(DE_rlog_data_2v3)

# Heatmap generation
if (nrow(DE_rlog_matrix_2v3) > 1) {
  heatmap.2(
    DE_rlog_matrix_2v3, 
    scale = "row",
    trace = "none",
    col = color_palette,
    main = "Heatmap entre le groupe 2 et le groupe 3",
    
    # Personalization 
    labRow = NA,
    cexCol = 1,        
    srtCol = 45,       
    margins = c(8, 10)
  )
} else {
  cat("Aucun gène significatif pour Groupe 2 vs Groupe 3\n")
}


#dev.off()
