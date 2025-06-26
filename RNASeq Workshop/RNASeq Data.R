

# Install required packages (if not already installed)
if (!require("DESeq2")) install.packages("DESeq2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggrepel")) install.packages("ggrepel")

# Load libraries
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)

PATH <- "c:/Users/mbzhab/Desktop/RNASeq Workshop"
setwd (PATH)

PATH <- "C:/Users/mbzhab/OneDrive - The University of Nottingham/R Workshops/RNASeq Workshop"
setwd (PATH)

# Read the 'Raw-counts.tsv'file
counts <- read.delim(paste(PATH, "/DATA/Raw-counts.tsv", sep=""))

# Read the 'Experiment-design.tsv' file
design <- read.delim(paste(PATH, "/DATA/Experiment-design.tsv", sep=""))

# Inspect data
dim (counts)
dim (design)

counts[1:5, 1:5]

head(design)

table (design$Tissue)
table (design$Tissue, design$Therapy)

# Subset of ExpDes to create metadata
metadata <- data.frame(Sample_ID = design$Run, Type = as.factor(design$Type))

head (metadata)

# Changing the rownames of the metadata to Sample IDs
row.names(metadata) <- metadata$Sample_ID
head (metadata)

# Chang the rownames of counts dataset to Ensemble ID
row.names(counts) <- counts$ENS_ID

# Remove the ENS_ID from the datafram
counts <- counts %>% select(-ENS_ID)
head (counts)

dim (counts)

# Check if there are any differences between values stored in two variables
identical (colnames(counts), rownames(metadata))
setdiff (colnames(counts), rownames(metadata))  
intersect (colnames(counts), rownames(metadata)) 

# Check if the Sample_IDs in both datasets are in the same orders
all (colnames(counts) == rownames(metadata))    

# Sort rownames of metadata and colnames of counts in ascending
metadata <- metadata[order(rownames(metadata)), ]
counts <- counts[, order(colnames(counts))]

all (colnames (counts) == rownames(metadata))    

#Perform DESeq analysis
dds <- DESeqDataSetFromMatrix (counts, metadata, ~Type)
dds <- DESeq (dds)

# Create DESeq results object using Benjamini-Hochberg correction
res <- results (object = dds, contrast = c ('Type', 'cancer', 'normal'),
                pAdjustMethod = 'BH', alpha = 0.05)

res
summary (res)

# Variance Stabilizing Transformation (VST)
vst <- vst (dds)

# Enhanced PCA Plot
pcaData <- plotPCA(vst, intgroup = "Type", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = Type,
    label = rownames(pcaData))) +
  geom_point(size = 4, shape = 19) +
  geom_text_repel(size = 3, max.overlaps = 10,
    show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic(base_size = 14) +
  ggtitle("PCA Plot of Samples")+
  theme(
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5), 
    panel.border = element_rect(color = "black",
      fill = NA, size = 0.8))

# Import annotation with Symbol names, gene position, and biotype
rel_ENS_to_SYM <- read.delim (paste (PATH, "/DATA/rel_ENS_SYM_chr.tsv", sep=""))
# Create a new dataframe containing the results
output <- as.data.frame (res)
head (output)

# Merge the results with the annotation dataframe
output <- merge(output, rel_ENS_to_SYM, by.x="row.names", by.y= "ENS_ID")
head (output)

# Sort the output on padj
output <- output [order (output$padj), ]
head (output)

# select protein coding genes and put them in output.coding
output.coding <- subset(output, Gene_type == "protein_coding")

# Identify significant Differentially Expressed Genes
signigficantDEGs <- subset (output.coding, padj < 0.05)
UP_regulated <- subset (signigficantDEGs, log2FoldChange > 0)
DOWN_regulated <- subset (signigficantDEGs, log2FoldChange < 0)

# Anotate genes in Volcano Plot 
output.coding$DiffExpressed <- "NO"
output.coding$DiffExpressed [output.coding$log2FoldChange > 0 & output.coding$padj < 0.0001] <- "UP"
output.coding$DiffExpressed [output.coding$log2FoldChange < 0 & output.coding$padj < 0.0001] <- "DOWN"

# Create a new column DElabel, which contains the name of DEGs
output.coding$DElabel <- NA
output.coding$DElabel [output.coding$SYMBOL %in% 
                         c(head(DOWN_regulated$SYMBOL, 5), head(UP_regulated$SYMBOL, 5))] <- 
  output.coding$SYMBOL[output.coding$SYMBOL %in%
                         c(head(DOWN_regulated$SYMBOL, 5), head(UP_regulated$SYMBOL, 5))]


# Draw the Volcano plot
ggplot(data = subset(output.coding, !is.na(padj)),
    aes(x = log2FoldChange, y = -log10(padj),
    color = DiffExpressed, label = DElabel)) +
  geom_point(alpha = 0.7) +
  theme_classic() +
  geom_text_repel(size = 3, max.overlaps = 10,
    show.legend = FALSE) +
  scale_color_manual(values = c("blue", "black", "red")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value", color = "Expression") +
  theme(
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5), 
    panel.border = element_rect(color = "black",
      fill = NA, size = 0.8))


# Density plot of log2 fold changes
ggplot(output.coding, aes(x = log2FoldChange,
    fill = DiffExpressed)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("blue", "red", "grey"),
    name = "Gene Expression") +
  theme_classic() +
  labs(title = "Density Plot of Log2 Fold Changes",
    x = "Log2 Fold Change", y = "Density", legend= "He") +
  theme(
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5), 
    panel.border = element_rect(color = "black",
     fill = NA, size = 0.8))
