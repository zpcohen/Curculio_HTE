## codon bias:
library(coRdon)

# 1. Load cDNA set
dnaCCall <- readSet(file="/home/zach/codon_bias/25_70_C13_CCB_PCN_Ribo.input")
CCall <- codonTable(dnaCCall)
genes <- getKO(CCall)[getlen(CCall) > 10]
summary(CCall)
head(names(dnaCCall))
setwd("/home/zach/codon_bias")

# 2. Get IDs
ribo_ids <- readLines("PCN_ribo_genes.list")
c1_hgt_ids <- readLines("Chr1_HGT.list")
c3_hgt_ids <- readLines("Chr3_HGT.list")
pcn_ids <- readLines("PCN_24k.list")
CCB_ids <- readLines("CCB_gene.list")
R_HGT <- readLines("70_fullheaders.list")
Rickettsia <- readLines("Ricketts_gene.list")

# 2. Initialize labels
labels <- rep("other", length(dnaCCall))
names(labels) <- names(dnaCCall)

table(labels)
# 3. Apply fuzzy matching to headers
label_by_ids <- function(names_vec, id_list) {
  vapply(names_vec, function(x) any(grepl(paste(id_list, collapse = "|"), x)), logical(1))
}

is_ribo <- names(dnaCCall) %in% ribo_ids
is_c1HGT <- names(dnaCCall) %in% c1_hgt_ids
is_c3HGT <- names(dnaCCall) %in% c3_hgt_ids
is_PCN <- names(dnaCCall) %in% pcn_ids
is_ccb <- names(dnaCCall) %in% CCB_ids
is_rHGT <- names(dnaCCall) %in% R_HGT
is_ricketts <- names(dnaCCall) %in% Rickettsia

length(labels)
# 3. Initialize a full label vector
labels <- rep("other", length(dnaCCall))
names(labels) <- names(dnaCCall)
table(labels)

# 4. Assign group labels by logical indexing
genes[is_ribo] <- "ribosomal"
genes[is_c1HGT] <- "PCN_Ccb"
genes[is_PCN] <- "PCN"
genes[is_c3HGT] <- "PCN_Ccb"
genes[is_ccb] <- "Ccb"

genes[is_rHGT] <- "PCN_RCa"
genes[is_ricketts] <- "RCa"

labels[is_ribo] <- "ribosomal"
labels[is_c1HGT] <- "PCN_Ccb"
labels[is_PCN] <- "PCN"
labels[is_c3HGT] <- "PCN_Ccb"
labels[is_ccb] <- "Ccb"
labels[is_rHGT] <- "PCN_Rca"
labels[is_ricketts] <- "Rca"

table(labels)
length(labels)
other_genes_idx <- which(labels == "other")


all(names(CCall) == names(labels))  # should return TRUE
length(labels)

milc_all <- MILC(CCall, subsets = list(
  ribosomal = labels == "ribosomal",
  C1_HGT = labels == "PCN_Ccb",
  C3_HGT = labels == "PCN_Ccb",
  CCB = labels == "Ccb",
  PCN = labels == "PCN",
  C_RCa =  labels == "PCN_Rca",
  RCa = labels == "Rca"
))

milc_df <- as.data.frame(milc_all)
head(milc_df)
milc_df$label <- labels
ggplot(milc_df, aes(x = PCN, y = ribosomal, color = label)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "Codon Bias in Ribosomal, HGT, PCN, and Other Genes",
    x = "MILC from Ribosomal CU",
    y = "MILC from Global CU"
  )

kruskal.test(ribosomal ~ label, data = milc_df)
# Kruskal-Wallis chi-squared = 927.74, df = 5, p-value < 2.2e-16

pairwise.wilcox.test(milc_df$ribosomal, milc_df$label, p.adjust.method = "BH")

#          Ccb     PCN_Ccb Rca    PCN_Rca PCN    
# PCN_Ccb   5.2e-11 -       -      -       -      
#  Rca       2.4e-07 0.0097  -      -       -      
#  PCN_Rca   < 2e-16 1.8e-13 0.0014 -       -      
#  PCN       < 2e-16 < 2e-16 0.0179 0.0732  -      
#  ribosomal 4.2e-10 0.1009  0.2975 2.7e-07 1.6e-14

milc_df$label <- factor(
  milc_df$label,
  levels = c("Ccb", "PCN_Ccb", "Rca", "PCN_Rca", "PCN", "ribosomal")
)
table(milc_df$label)
unique(milc_df$label)
#violin plot
ggplot(milc_df, aes(x = label, y = ribosomal, fill = label)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, position = position_dodge(width = 0.9)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "MILC vs Ribosomal Codon Usage",
    x = "Gene Group",
    y = "MILC (relative to ribosomal CU)"
  )

# PCA with centroid
codon_mat <- codonCounts(CCall)

###
gene_ids <- (CCall@ID)
head(gene_ids)
length(gene_ids)
rownames(codon_mat) <- gene_ids
length(rownames(codon_mat))
head(rownames(codon_mat))

codon_freqs <- sweep(codon_mat, 1, rowSums(codon_mat), FUN = "/")
length(codon_freqs)
codon_freqs <- as.data.frame(codon_freqs)

codon_freqs$label <- labels
dim(codon_freqs)
head(codon_freqs$label)
pca <- prcomp(codon_freqs[, -ncol(codon_freqs)], scale. = TRUE)

pca_df <- data.frame(pca$x[, 1:2])  # Keep only PC1 and PC2 for plotting
pca_df$label <- codon_freqs$label  # Add gene labels

# Cc = "#4E79A7",
#Cc_Ccb = "#E15759",
#Cc_R = "#59A14F")) +

label_colors <- c(
  "Ccb"       = "#F8766D",
  "PCN Ccb"    = "#E15759",
  "Rca"    = "#00BA38",
  "PCN_Rca"       = "#59A14F",
  "PCN" = "#4E79A7",
  "ribosomal" = "#F781BF"
)

centroids <- pca_df %>%
  group_by(label) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

ggplot(pca_df, aes(x = PC1, y = PC2, color = label)) +
  geom_point(data = subset(pca_df, label == "PCN"),
             aes(x = PC1, y = PC2, color = label),
             alpha = 0.001, size = 1) + 
  geom_point(data = subset(pca_df, label != "PCN"),
             aes(x = PC1, y = PC2, color = label),
             alpha = 0.7, size = 1.5) +
  geom_point(alpha = 0.6) +
  geom_point(data = centroids, aes(x = PC1, y = PC2, color = label),
             shape = 21, fill = "white", size = 4, stroke = 1.2) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(), aspect.ratio = 1) +
  labs(
    title = "PCA of Codon Usage Frequencies",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "% variance)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "% variance)")
  ) +
  scale_color_manual(values = label_colors) 


label_colors <- c(
  "Ccb" = "#F8766D",        # red
  "PCN" = "#4E79A7",        # blue
  "PCN_Ccb" = "#E7872B",    # gold
  "PCN_Rca" = "#59A14F",    # green
  "Rca" = "#00C1AA",        # teal
  "ribosomal" = "#FF61CC"   # pink
)

ggplot() +
  # Transparent PCN group
  geom_point(data = subset(pca_df, label == "PCN"),
             aes(x = PC1, y = PC2, color = label),
             alpha = 0.3, size = 1) +  # lower alpha for PCN
  
  geom_point(data = subset(pca_df, label == "Ccb"),
             aes(x = PC1, y = PC2, color = label),
             alpha = 0.6, size = 1) +  # lower alpha for Ccb
  # Other groups with default alpha
  geom_point(data = subset(pca_df, !(label %in% c("PCN", "Ccb"))),
             aes(x = PC1, y = PC2, color = label),
             alpha = 0.7, size = 1.5) +
  
  # Add centroids if needed
  geom_point(data = centroids, aes(x = PC1, y = PC2, fill = label),
             size = 5, shape = 21, color = "white", stroke = 1.5) +
  
  scale_color_manual(values = palette) +  # your custom color palette
  scale_fill_manual(values = palette) +
  labs(
    title = "PCA of Codon Usage Frequencies",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "% variance)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "% variance)")
  ) +
  theme_minimal(base_size = 14)
library(dplyr)

centroids <- pca_df %>%
  group_by(label) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

ggplot(pca_df, aes(x = PC1, y = PC2, color = label)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_point(data = centroids, aes(x = PC1, y = PC2, color = label),
             shape = 21, fill = "white", size = 4, stroke = 1.2) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(), aspect.ratio = 1) +
  labs(title = "PCA of Codon Usage Frequencies (with group centroids)") +
  scale_color_manual(values = label_colors)


## intron length R code:
df_all <- read.table("/home/zach/codon_bias/All_input_Feb2626.l2", header = T)


sample_counts <- df_all %>%
  group_by(group) %>%
  summarise(n = n(), .groups = "drop")

p <- ggplot(df_all, aes(x = group, y = intron_length, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values =c(
    Cc = "#4E79A7",
    Cc_Ccb = "#E15759",
    Cc_R = "#59A14F")) +
    scale_y_continuous(trans = "log10") +
 # geom_text(
#    data = sample_counts,
#    aes(x = group, y = 1e7, label = paste0("n = ", n)),  # adjust y = position as needed
#    inherit.aes = FALSE,
#    size = 6
#  ) +
  theme_minimal(base_size = 16) +
  theme(panel.grid = element_blank())+
  theme(aspect.ratio = 1) +
  labs(
    title = "Intron Length Distribution by Gene Group",
    x = "Group",
    y = "Total Intron Length per Gene"
  ) 
p + stat_compare_means(
  comparisons = list(
    c("Cc", "Cc_Ccb"),
    c("Cc", "Cc_R"),
    c("Cc_Ccb", "Cc_R")
  ),
  method = "wilcox.test",
  label = "p.signif",  # or "p.format" for numeric values
  size = 6
  )

wilcox.test(intron_length ~ group, data = df_all)

