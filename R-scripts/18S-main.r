###################### Libraries ######################
library(tidyverse) 
library(qiime2R) # Integrating QIIME2 and R for data visualization and analysis using qiime2R.
library(phyloseq) # Handling and analysis of high-throughput microbiome census data.
library(ggh4x) # ggh4x package is a ggplot2 extension package.
library(cowplot) # Simple add-on to ggplot2; provides additional features to improve graphic quality.  
library(ggpubr) # Provides some easy-to-use functions to customize ggplot2 graphics.
library(decontam) # Decontam package is used to identify and remove contaminants from microbiome data.
library(fantaxtic) # Provides functions to work with taxonomic data in R.
library(microbiome) # Tools for microbiome analysis
library(gridExtra) # Provides functions to arrange multiple grid-based figures on a single page.
library(microViz) # Provides functions for microbiome data visualization and analysis.
library(vegan) # Community ecology package; it can be use for diverse multivariate analysis.
library(dplyr) # A tool for working/manipulating dataframes

##################### Importing and cleaning data ####
# Imported taxonomy, dada2 table, rooted tree, and metadata from local storage
taxonomy_18s <- qiime2R::read_qza("18s/QIIME2 Outputs/18S-rep-sequences-taxonomy.qza") 

table_18s <- qiime2R::read_qza("18s/QIIME2 Outputs/18S-dada2-table.qza")

tree_18s <- qiime2R::read_qza("18s/QIIME2 Outputs/rooted-18S-tree.qza") 

metadata_18s <- read.csv("Raw Data/seagrasses_metadata_habitat_location.csv",
                     sep = ",", header = T, row.names = 1, quote = "")

refseq <- qiime2R::read_qza("18s/QIIME2 Outputs/denoised-rep-sequences.qza") 

refseq_df <- refseq$data

physeq_seq <- refseq(refseq_df)

# Extracted taxonomy data from taxonomy file
taxonomy_info_18s <- taxonomy_18s$data

# Parsed the taxonomy into 22 different taxonomic levels
# V1 is Domain, V2 is Kingdom, V14 is Phylum, V15 is Class, V17 is Order, V21 is Family, and V22 is Genus
parse_taxonomy_18s <- taxonomy_info_18s %>% separate(Taxon, sep =";", into = c("V1","V2","V3","V4","V5","V6","V7","V8",
                                                                                       "V9","V10","V11","V12","V13","V14","V15","V16","V17",
                                                                                       "V18","V19","V20","V21","V22")) 
# Made Feature IDs the row names and removed the Feature.ID column 
rownames(parse_taxonomy_18s) <- parse_taxonomy_18s$Feature.ID 

parse_taxonomy_18s$Feature.ID <- NULL 

# Removed any digits or characters from the beginning of taxon names
# Empty cells were classified as NA
# Replaced all NAs with "Unassigned"
parse_taxonomy_18s[] <- lapply(parse_taxonomy_18s, function(x) gsub("D_\\d+__", "", x))

parse_taxonomy_18s[parse_taxonomy_18s==""] <- NA

parse_taxonomy_18s <- parse_taxonomy_18s %>% replace(is.na(.), "Unassigned")

write.csv(parse_taxonomy_18s,"18s/Clean Data/taxonomy.csv")

# Transformed parsed taxonomy, table, and rooted tree into phyloseq objects
TAX_18s <- phyloseq::tax_table(as.matrix(parse_taxonomy_18s))

OTUMAT_18s <- table_18s$data

OTU_18s <- otu_table(OTUMAT_18s, taxa_are_rows = TRUE, replace_na(0))

otu_18s_export <- as.data.frame(OTU_18s)

write.csv(otu_18s_export,"18s_asv_table.csv")

TREE_18s <- tree_18s$data

# Created phyloseq and removed any samples that had less than 1 read
phylo_18s <- phyloseq(OTU_18s, sample_data(metadata_18s), TAX_18s, TREE_18s) # 5470 taxa and 196 samples

phylo_18s <- merge_phyloseq(phylo_18s, physeq_seq)

phylo_18s <- phylo_18s %>% 
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .) # 5470 taxa and 193 samples

# Removal of susupicious samples and ordination 
to_remove4 <- c("18S-bodega-bay_CC.B.1.1.RS") 

to_remove5 <- c("18S-bodega-bay_CC.B.3.3.RS") 

phylo_18s <- prune_samples(!(sample_names(phylo_18s) %in% to_remove4), phylo_18s) # 5470 taxa and 192 samples

phylo_18s <- prune_samples(!(sample_names(phylo_18s) %in% to_remove5), phylo_18s) # 5470 taxa and 191 samples

##################### Using Decontam for filering out contaminants ####
sample_data(phylo_18s)$is.neg <- sample_data(phylo_18s)$Habitat == "NegtCtrl" | sample_data(phylo_18s)$Habitat == "Mock" | sample_data(phylo_18s)$Habitat == "Blank" # create a sample-variable for contaminants
phylo_contaminants_18s <- isContaminant(phylo_18s, method = "prevalence", neg="is.neg", threshold=0.5, detailed = TRUE, normalize = TRUE) # detect contaminants based on control samples and their ASV prevalance
table(phylo_contaminants_18s$contaminant) # check number of ASVs that are contaminents (41)

# Make phyloseq object of presence-absence in negative controls and true samples
phylo_contaminants.pa_18s <- transform_sample_counts(phylo_18s, function(abund) 1 * (abund > 0)) # convert phyloseq table to presence-absence

ps.pa.neg_18s <- prune_samples(sample_data(phylo_contaminants.pa_18s)$Habitat == "Control" | sample_data(phylo_contaminants.pa_18s)$Habitat == "Blank" | 
                                 sample_data(phylo_contaminants.pa_18s)$Habitat == "Mock",phylo_contaminants.pa_18s) # identify controls

ps.pa.pos_18s <- prune_samples(sample_data(phylo_contaminants.pa_18s)$Habitat != "Control" | sample_data(phylo_contaminants.pa_18s)$Habitat != "Blank" |
                                 sample_data(phylo_contaminants.pa_18s)$Habitat != "Mock", phylo_contaminants.pa_18s) # identify samples

df.pa_18s <- data.frame(pa.pos=taxa_sums(ps.pa.pos_18s), pa.neg=taxa_sums(ps.pa.neg_18s), contaminant=phylo_contaminants_18s$contaminant) # convert into a dataframe

# Make phyloseq object of presence-absence in negative controls and true samples
ggplot(data=df.pa_18s, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

phylo_18s <- prune_taxa(!phylo_contaminants_18s$contaminant, phylo_18s) # remove ASVs identified as decontaminants from the dataset (5429 taxa and 191 samples)

df_phylo_18s <- psmelt(phylo_18s)

clan_18s_export <- as.data.frame(otu_table(phylo_18s))

write.csv(clan_18s_export, "clan_18s_export.csv")

###################### Making edits to original phyloseq ####################### 
# Removed the group Charophyta (sea grasses)
phylo_18s <- subset_taxa(phylo_18s, V4 != "Charophyta") # 5294 taxa and 191 samples

phylo_18s <- subset_taxa(phylo_18s, V1 != "Unassigned") # 4165 taxa and 191 samples

phylo_18s <- subset_taxa(phylo_18s, !grepl("Ulva", V6, ignore.case = TRUE)) # 4150 taxa and 191 samples

# Subset phyloseq by Habitat
# This removed any mocks, blanks, or controls that were left in the original phyloseq
# Normalized the phyloseq for relative abundance 
phylo_18s <- phylo_18s %>% subset_samples(Site %in% c("Mason's Marina", "Campbell Cove", "Westside Park")) # 4150 taxa and 165 samples

phylo_18s <- phylo_18s %>% 
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)

phylo_normalized_18s <- microbiome::transform(phylo_18s, "compositional")

# Subset the samples that were Raw Sediment sample type
# Normalized the Raw Sediment samples for relative abundance 
phylo_raw_18s <- phylo_18s %>% subset_samples(SampleType %in% c("RawSediment")) # 4150 taxa and 81 samples

phylo_normalized_raw_18s <- microbiome::transform(phylo_raw_18s, "compositional")

# Subset the samples that were Ludox sample type
# Normalized the Ludox samples for relative abundance 
phylo_ludox_18s <- phylo_18s %>% subset_samples(SampleType %in% c("Ludox")) # 4150 taxa and 84 samples

phylo_normalized_ludox_18s <- microbiome::transform(phylo_ludox_18s, "compositional")

# Subset the samples for only Nematodes from Ludox samples
phylo_nematoda_18s <- subset_taxa(phylo_ludox_18s, V14 %in% c("Nematoda")) # 772 taxa and 84 samples

phylo_normalized_nematoda_18s <- microbiome::transform(phylo_nematoda_18s, "compositional")

df_phylo_ludox <- psmelt(phylo_ludox_18s)

df_phylo_raw <- psmelt(phylo_raw_18s)

# Site specific phyloseq objects
CC_phylo_18s <- phylo_18s %>% subset_samples(Site %in% ("Campbell Cove"))
CC_phylo_18s_normalized <- microbiome::transform(CC_phylo_18s, "compositional")

WP_phylo_18s <- phylo_18s %>% subset_samples(Site %in% ("Westside Park"))
WP_phylo_18s_normalized <- microbiome::transform(WP_phylo_18s, "compositional")

MM_phylo_18s <- phylo_18s %>% subset_samples(Site %in% ("Mason's Marina"))
MM_phylo_18s_normalized <- microbiome::transform(MM_phylo_18s, "compositional")

# Site specific phyloseq objects by Sample Type
CC_phylo_18s_raw <- phylo_raw_18s %>% subset_samples(Site %in% ("Campbell Cove"))
CC_phylo_18s_normalized_raw <- microbiome::transform(CC_phylo_18s_raw, "compositional")

WP_phylo_18s_raw <- phylo_raw_18s %>% subset_samples(Site %in% ("Westside Park"))
WP_phylo_18s_normalized_raw <- microbiome::transform(WP_phylo_18s_raw, "compositional")

MM_phylo_18s_raw <- phylo_raw_18s %>% subset_samples(Site %in% ("Mason's Marina"))
MM_phylo_18s_normalized_raw <- microbiome::transform(MM_phylo_18s_raw, "compositional")

CC_phylo_18s_ludox <- phylo_ludox_18s %>% subset_samples(Site %in% ("Campbell Cove"))
CC_phylo_18s_normalized_ludox <- microbiome::transform(CC_phylo_18s_ludox, "compositional")

WP_phylo_18s_ludox <- phylo_ludox_18s %>% subset_samples(Site %in% ("Westside Park"))
WP_phylo_18s_normalized_ludox <- microbiome::transform(WP_phylo_18s_ludox, "compositional")

MM_phylo_18s_ludox <- phylo_ludox_18s %>% subset_samples(Site %in% ("Mason's Marina"))
MM_phylo_18s_normalized_ludox <- microbiome::transform(MM_phylo_18s_ludox, "compositional")

# Site specific phyloseq objects for Nematodes only
CC_phylo_18s_nematode <- phylo_nematoda_18s %>% subset_samples(Site %in% ("Campbell Cove"))
CC_phylo_18s_normalized_nematode <- microbiome::transform(CC_phylo_18s_nematode, "compositional")

WP_phylo_18s_nematode <- phylo_nematoda_18s %>% subset_samples(Site %in% ("Westside Park"))
WP_phylo_18s_normalized_nematode <- microbiome::transform(WP_phylo_18s_nematode, "compositional")

MM_phylo_18s_nematode <- phylo_nematoda_18s %>% subset_samples(Site %in% ("Mason's Marina"))
MM_phylo_18s_normalized_nematode <- microbiome::transform(MM_phylo_18s_nematode, "compositional")

phylo_18s_all_samples <- microbiome::transform(phylo_ludox_18s, "compositional")
simple_ord_18s <- ordinate(phylo_18s_all_samples, "NMDS", "bray")
plot_ordination(phylo_18s_all_samples, simple_ord_18s, type="samples", color="Site", shape="Habitat")#, label = "Description")