# Libraries
library(phyloseq)
library(pheatmap)
library(Hmisc)
library(corrplot)
library(microViz)
library(igraph)


########### Making Correlation matrix ##################
# Make phylo objects 
bb_18s <- phyloseq(otu_table(phylo_nematoda_18s),tax_table(phylo_nematoda_18s), sample_data(phylo_nematoda_18s)) # 772 taxa and 84 samples

bb_16s <- phyloseq(otu_table(phylo_16s),tax_table(phylo_16s), sample_data(phylo_16s)) # 9013 taxa and 82 samples

# Remove samples that are not found in 16s
bb_18s <- subset_samples(bb_18s, Description != "CC.SG.3.1.LUD.Ludox")

bb_18s <- subset_samples(bb_18s, Description != "MM.SG.1.3.LUD.Ludox") # 772 taxa and 82 samples

# Edit sample names of nematode metadata to match 16s metadata
variable_cor_18s <- as.data.frame(sample_data(bb_18s))

variable_cor_18s$samples <- row.names(variable_cor_18s) 

variable_cor_18s[] <- lapply(variable_cor_18s, function(x) gsub(".LUD", "", x))

variable_cor_18s[] <- lapply(variable_cor_18s, function(x) gsub("18S-bodega-bay", "16S-bodega-bay-fecal-experiment", x))

variable_cor_18s <- variable_cor_18s[, c( 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)]

rownames(variable_cor_18s) <- variable_cor_18s$samples 

variable_cor_18s$samples <- NULL 

sample_data(bb_18s) <- variable_cor_18s


# Edit samples names of nematode otu table to match 16s dataset
otus_cor_18s <- otu_table(bb_18s)

write.csv(otus_cor_18s,"otus_cor_18s.csv")

otus_cor_18s <- read.csv("otus_cor_18s.csv", sep = ",", header = T, row.names = 1, quote = "")

colnames(otus_cor_18s) <- gsub("X.", "", colnames(otus_cor_18s))

colnames(otus_cor_18s) <- gsub("18S.bodega.bay", "16S-bodega-bay-fecal-experiment", colnames(otus_cor_18s))

colnames(otus_cor_18s) <- gsub(".LUD.", "", colnames(otus_cor_18s))

otus_cor_18s <- otu_table(otus_cor_18s, taxa_are_rows = TRUE)

# Generate new phyloseq object for nematodes
bb_18s_2 <- phyloseq(otus_cor_18s,tax_table(phylo_nematoda_18s), variable_cor_18s)

# Filter for 5% and 10% prevalence level
bb_18s_05 <- tax_filter(bb_18s_2, min_prevalence = .05) # 117 taxa and 82 samples

bb_16s_05 <- tax_filter(bb_16s, min_prevalence = .05) # 1806 taxa and 82 samples

bb_18s_10 <- tax_filter(bb_18s_2, min_prevalence = .1) # 68 taxa and 82 samples

bb_16s_10 <- tax_filter(bb_16s, min_prevalence = .1) # 1073 taxa and 82 samples

bb_18s_15 <- tax_filter(bb_18s_2, min_prevalence = .15) # 44 taxa and 82 samples

bb_16s_15 <- tax_filter(bb_16s, min_prevalence = .15) # 769 taxa and 82 samples

bb_18s_20 <- tax_filter(bb_18s_2, min_prevalence = .20) # 433 taxa and 82 samples

bb_16s_20 <- tax_filter(bb_16s, min_prevalence = .20) # 585 taxa and 82 samples

bb_18s_25 <- tax_filter(bb_18s_2, min_prevalence = .25) # 23 taxa and 82 samples

bb_16s_25 <- tax_filter(bb_16s, min_prevalence = .25) # 4376 taxa and 82 samples



# Normalize
bb_18s_2_normalized <- microbiome::transform(bb_18s_05, "compositional") 

bb_16s_normalized <- microbiome::transform(bb_16s_05, "compositional") 

# Fix unknowns
bb_18s_2_normalized <- tax_fix(bb_18s_2_normalized, unknowns = "Unassigned")

bb_16s_2_normalized <- tax_fix(bb_16s_normalized, unknowns = c("Unassigned", "uncultured", "uncultured bacterium", "uncultured organism",
                                                               "uncultured marine bacterium", "uncultured archaeon"))

# Agglomerate taxa at genus level 
bb_18s_2_glom <- tax_glom(bb_18s_2_normalized, taxrank = "V22") # 41 taxa

bb_16s_2_glom <- tax_glom(bb_16s_2_normalized, taxrank = "V6") # 199 taxa


# Extract OTU tables (taxa as columns, samples as rows)
otu1 <- as(otu_table(bb_18s_2_glom), "matrix") %>% t() %>% as.data.frame()
otu2 <- as(otu_table(bb_16s_2_glom), "matrix") %>% t() %>% as.data.frame()

# Convert OTU tables to matrices (samples × taxa)
otu1_matrix <- as.matrix(otu1)  # from ps1 (samples × taxa)
otu2_matrix <- as.matrix(otu2)  # from ps2 (samples × taxa)

# Compute correlations and p-values
cor_test <- rcorr(otu1_matrix, otu2_matrix, type = "spearman")  # or "pearson"

# Extract correlation coefficients and p-values
cor_matrix_coe <- cor_test$r[1:ncol(otu1), (ncol(otu1)+1):ncol(cor_test$r)]
p_matrix <- cor_test$P[1:ncol(otu1), (ncol(otu1)+1):ncol(cor_test$r)]

# correct *p*-values to control false discoveries (e.g., Benjamini-Hochberg FDR)
p_adjusted <- matrix(p.adjust(p_matrix, method = "fdr"), nrow = nrow(p_matrix))
rownames(p_adjusted) <- rownames(cor_matrix_coe)
colnames(p_adjusted) <- colnames(cor_matrix_coe)

# Keep significant
significant_cor <- cor_matrix_coe
significant_cor[p_adjusted > 0.05] <- NA  # Mask non-significant correlations

threshold <- 0.4
significant_cor <- significant_cor
significant_cor[abs(significant_cor) < threshold] <- NA


# Extract taxonomy tables
tax1 <- as.data.frame(tax_table(bb_18s_2_glom))  
tax2 <- as.data.frame(tax_table(bb_16s_2_glom))  

# Example: Create a label like "Genus_species" for each OTU
tax1$label <- paste(tax1$V22)  
tax2$label <- paste(tax2$V6)  

# If Genus is NA, use higher rank (e.g., Family)
tax1$label <- ifelse(is.na(tax1$V22), tax1$V20, tax1$label)  
tax2$label <- ifelse(is.na(tax2$V6), tax2$V5, tax2$label)  

# Assuming cor_matrix was created from otu1 (rows) and otu2 (columns)
rownames(significant_cor) <- tax1[rownames(significant_cor), "label"]  
colnames(significant_cor) <- tax2[colnames(significant_cor), "label"]  

# Transpose 
significant_cor_t <- t(significant_cor)

# Remove rows where all elements are zero
mat <- significant_cor[rowSums(!is.na(significant_cor)) > 0, ]

# Remove columns where all elements are zero
mat <- mat[, colSums(!is.na(mat)) > 0]

# Plot correlations that are significant
corrplot(mat, is.corr = TRUE, method = 'color', addCoef.col = 'black', addgrid.col = 'grey',
         tl.srt = 45, tl.col = "black", col = COL2('BrBG'),
         na.label = 'square', na.label.col = "white", number.cex = 0.2, tl.cex = 0.3)

############ Exporting as igraph ###############

matrix_for_export <- mat

matrix_for_export[is.na(matrix_for_export)] <- 0

matrix_export_df <- as.data.frame(matrix_for_export)

write.csv(matrix_export_df, "co-occurence_matrix_04_threshold_01_sig.csv")

bip_graph <- graph_from_incidence_matrix(
  matrix_for_export,
  weighted = TRUE,
  mode = "all",
)

plot(bip_graph, layout = layout_on_sphere)

# Edge list
edges <- as_data_frame(bip_graph, what = "edges")
write.csv(edges, "edges.csv", row.names = FALSE)

# Node list
nodes <- as_data_frame(bip_graph, what = "vertices")
write.csv(nodes, "nodes.csv", row.names = FALSE)


# Project to one-mode network (e.g., document similarities)
proj <- bipartite_projection(bip_graph)
doc_network <- proj$proj1

# Export to Gephi (GraphML format)
write_graph(bip_graph, "all_at_04_threshold.graphml", format = "graphml")

write_graph(bip_graph, "new_net.graphml", format = "graphml")
