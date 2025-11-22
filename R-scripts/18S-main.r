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

###################### Importing and cleaning data ######################
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

###################### Using Decontam for filering out contaminants ######################
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

###################### Making edits to original phyloseq ###################### 
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

###################### Plotting for Relative Abundance Top 20 and Top 10 Raw Sediment at Phylum ######################
# Raw Sediment
# Function to collapse a certain number of taxa into category others
merge_top20_18s_raw <- function(dataframe_phylo, top=19){
  transformed <- transform_sample_counts(dataframe_phylo, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:23 if there are species level
  }
  return(merged)
}

# Agglomerated taxa down to V14 rank 
# Use tax_fix function to assign taxa names to lower ranks that are unknown
temp_phylo_18s_raw <- phylo_raw_18s

fix_phylo_18s_raw <- tax_fix(temp_phylo_18s_raw, unknowns = c("Unassigned", "uncultured", "Unkown", "uncultured eukaryote", "Unknown Family",
                                                              "uncultured Cercozoa", "uncultured dinoflagellate"))

glom_18s_raw <- tax_glom(fix_phylo_18s_raw, taxrank = "V14")

# Run function on phyloseq object and make into data frame
phy_18_raw_top20_V14 <- merge_top20_18s_raw(glom_18s_raw, top=19)

phy_18_raw_top20_V14_df <- psmelt(phy_18_raw_top20_V14)

write.csv(phy_18_raw_top20_V14_df, "phy_18_raw_top20_V14_df.csv")

# Add common factors to use for plotting
phy_18_raw_top20_agr = aggregate(Abundance~Sample+Site+Habitat+V14, data=phy_18_raw_top20_V14_df, FUN=mean) 
unique(phy_18_raw_top20_agr$V14)

# Organize color scales and factors to plot top10 taxa as barplots
# List of 20 distinct colors
colors_top20 <- c("#A6CEE3", "#579CC7", "#3688AD", "#8BC395", "#89CB6C", "#40A635", "#919D5F", "#F99392", "#EB494A","#F79C5D",
                  "#FDA746", "#FE8205", "#E39970", "#BFA5CF", "#8861AC", "#917099", "#E7E099", "#DEB969", "#B15928", "gray")

colors_top20_2 <- c("#9a6324", "#46f0f0", "#1F78B4", "#aaffc3", "#e6beff",  "#33A02C", "#4363d8", "#008080", "#FB9A99", "#e6194b",
                  "#f032e6", "#bcf60c", "#fabebe","#f58231", "#ffe119", "#6A3D9A", "#000075", "#fffac8", "#800000","gray")


# Create new labels for important factors
habitat.labs <- c("Bare Sediment", "Sea Grass")
names(habitat.labs) <- c("Bare Sediment", "Sea Grass")

site.labs <- c("Campbell Cove", "Westside Park", "Mason's Marina")
names(site.labs) <- c("Campbell Cove", "Westside Park", "Mason's Marina")

# Put "Others" to the final of the Phylum list - top 20
unique(phy_18_raw_top20_agr$V14)
phy_18_raw_top20_agr$V14 <- factor(phy_18_raw_top20_agr$V14,
                              levels = c("Ansanella natalensis V9", "Biecheleria V8", "Brachiopoda","Bysmatrum subsalsum V8", "Cocconeis placentula V9",
                                         "Crustacea", "eukaryote marine clone ME1-22 V6", "Gonyaulax spinifera V9", "Haslea ostrearia V9", "Lankesteria V8",
                                         "Lecudina V8", "Maullinia V6", "Minidiscus sp. V9","Nematoda", "Opephora sp. s0357 V8", 
                                         "Polychaeta", "Rhabditophora","uncultured Cercozoa V7", "uncultured dinoflagellate V9", "Others"))

# Reorder Site levels
phy_18_raw_top20_agr$Site = factor(phy_18_raw_top20_agr$Site, levels=c("Campbell Cove","Westside Park","Mason's Marina"))

# Plot by site - Phylum level top 20
taxonomy_bar_18s_raw_top20_phylum <- ggplot(phy_18_raw_top20_agr, aes(x = Sample, y = Abundance, fill = V14)) +
  facet_nested(. ~ Site+Habitat, scales = "free",
               labeller = labeller(Site = site.labs, Habitat = habitat.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20, name = "Raw Sediment Phylum") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, vjust = 0.5, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12)) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) + # Format facet grid title 
  scale_x_discrete(label = function(x) stringr::str_replace(x, "18S-bodega-bay_","")) 

taxonomy_bar_18s_raw_top20_phylum

# Function for top 10
merge_top10_18s_raw <- function(phylo_raw_18s, top=9){
  transformed <- transform_sample_counts(phylo_raw_18s, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:23 if there are species level
  }
  return(merged)
}


# Run function on phyloseq object and make into data frame
phy_18_raw_top10_V14 <- merge_top20_18s_raw(glom_18s_raw, top=9)

phy_18_raw_top10_V14_df <- psmelt(phy_18_raw_top10_V14)

# Add common factors to use for plotting
phy_18_raw_top10_agr = aggregate(Abundance~Sample+Site+Habitat+V14, data=phy_18_raw_top10_V14_df, FUN=mean) 
unique(phy_18_raw_top10_agr$V14)

# Organize color scales and factors to plot top10 taxa as barplots
# List of 10 distinct colors
colors_top10_2 <- c("#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#f5f5f5", "#c7eae5", "#80cdc1", "#35978f", "#01665e", "gray")

# Put "Others" to the final of the Phylum list - top 10
unique(phy_18_raw_top10_agr$V14)
phy_18_raw_top10_agr$V14 <- factor(phy_18_raw_top10_agr$V14,
                                   levels = c("Ansanella natalensis V9","Crustacea", "eukaryote marine clone ME1-22 V6", "Gonyaulax spinifera V9","Lankesteria V8",
                                              "Maullinia V6", "Minidiscus sp. V9","Nematoda",
                                              "Polychaeta","Others"))

# Reorder Site levels
phy_18_raw_top10_agr$Site = factor(phy_18_raw_top10_agr$Site, levels=c("Campbell Cove","Westside Park","Mason's Marina"))

# Plot by site - Class level top 10
taxonomy_bar_18s_raw_top10_phylum <- ggplot(phy_18_raw_top10_agr, aes(x = Sample, y = Abundance, fill = V14)) +
  facet_nested(. ~ Site+Habitat, scales = "free",
               labeller = labeller(Site = site.labs, Habitat = habitat.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top10_2, name = "Raw Sediment Phylum") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, vjust = 0.5, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12)) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) + # Format facet grid title 
  scale_x_discrete(label = function(x) stringr::str_replace(x, "18S-bodega-bay_","")) 

taxonomy_bar_18s_raw_top10_phylum

##################### Plotting for Relative Abundance Top 20 and Top 10 Ludox at Phylum ######################
# Function to collapse a certain number of taxa into category others
merge_top20_18s_ludox <- function(phylo_ludox_18s, top=19){
  transformed <- transform_sample_counts(phylo_ludox_18s, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:23 if there are species level
  }
  return(merged)
}

# Agglomerated taxa down to V14 rank 
# Use tax_fix function to assign taxa names to lower ranks that are unknown
temp_phylo_18s_ludox <- phylo_ludox_18s

fix_phylo_18s_ludox <- tax_fix(temp_phylo_18s_ludox, unknowns = c("Unassigned", "uncultured", "Unkown", "uncultured eukaryote", "Unknown Family"))

glom_18s_ludox <- tax_glom(fix_phylo_18s_ludox, taxrank = "V14")

# Run function on phyloseq object and make into data frame
phy_18_ludox_top20_V14 <- merge_top20_18s_ludox(glom_18s_ludox, top=19)

phy_18_ludox_top20_V14_df <- psmelt(phy_18_ludox_top20_V14)

# Add common factors to use for plotting
phy_18_ludox_top20_agr = aggregate(Abundance~Sample+Site+Habitat+V14, data=phy_18_ludox_top20_V14_df, FUN=mean) 
unique(phy_18_ludox_top20_agr$V14)

# Put "Others" to the final of the Phylum list - top 20
phy_18_ludox_top20_agr$V14 <- factor(phy_18_ludox_top20_agr$V14,
                                     levels = c("Alexandrium sp. CCMP1911 V9", "Alexandrium V8", "Brachiopoda", "Crustacea", "Gloiopeltis furcata V6",
                                                "Gonyaulax spinifera V9", "Gymnodinium sp. MUCC284 V8", "Halodaphnea panulirata V6", "Lankesteria V8",
                                                "Lecudina phyllochaetopteri V9", "Maullinia V6", "Oligotrichia sp. EP-2016a V9", "Polychaeta","Nematoda",
                                                "Protoperidinium conicum V9", "Protoperidinium leonis V9", "Rhabditophora", "Salispina spinosa var. spinosa V6",
                                                "Tribonema marinum V8", "Others"))
                                     
# Reorder Site levels
phy_18_ludox_top20_agr$Site = factor(phy_18_ludox_top20_agr$Site, levels=c("Campbell Cove","Westside Park","Mason's Marina"))

# Plot by site - Class level top 10
  taxonomy_bar_18s_ludox_top20_phylum <- ggplot(, aes(x = Sample, y = Abundance, fill = V14)) +
  facet_nested(. ~ Site+Habitat, scales = "free",
               labeller = labeller(Site = site.labs, Habitat = habitat.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20, name = "Ludox Phylum") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, vjust = 0.5, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12)) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) + # Format facet grid title 
  scale_x_discrete(label = function(x) stringr::str_replace(x, "18S-bodega-bay_","")) 


taxonomy_bar_18s_ludox_top20_phylum





merge_top10_18s_ludox <- function(phylo_ludox_18s, top=9){
  transformed <- transform_sample_counts(phylo_ludox_18s, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:23 if there are species level
  }
  return(merged)
}


# Run function on phyloseq object and make into data frame
phy_18_ludox_top10_V14 <- merge_top20_18s_ludox(glom_18s_ludox, top=9)

phy_18_ludox_top10_V14_df <- psmelt(phy_18_ludox_top10_V14)

# Add common factors to use for plotting
phy_18_ludox_top10_agr = aggregate(Abundance~Sample+Site+Habitat+V14, data=phy_18_ludox_top10_V14_df, FUN=mean) 
unique(phy_18_ludox_top10_agr$V14)

# Put "Others" to the final of the Phylum list - top 10
phy_18_ludox_top10_agr$V14 <- factor(phy_18_ludox_top10_agr$V14,
                                     levels = c("Crustacea", "Gonyaulax spinifera V9","Lecudina phyllochaetopteri V9",
                                                "Maullinia V6","Polychaeta", "Protoperidinium conicum V9", "Protoperidinium leonis V9", "Nematoda",
                                                "Rhabditophora","Others"))

# Reorder Site levels
phy_18_ludox_top10_agr$Site = factor(phy_18_ludox_top10_agr$Site, levels=c("Campbell Cove","Westside Park","Mason's Marina"))

# Plot by site - Class level top 10
taxonomy_bar_18s_ludox_top10_phylum <- ggplot(phy_18_ludox_top10_agr, aes(x = Sample, y = Abundance, fill = V14)) +
  facet_nested(. ~ Site+Habitat, scales = "free",
               labeller = labeller(Site = site.labs, Habitat = habitat.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top10_2, name = "Ludox Phylum") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, vjust = 0.5, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12)) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) + # Format facet grid title 
  scale_x_discrete(label = function(x) stringr::str_replace(x, "18S-bodega-bay_","")) 

taxonomy_bar_18s_ludox_top10_phylum
