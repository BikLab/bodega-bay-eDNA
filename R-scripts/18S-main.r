###################### Libraries ######################
library(tidyverse) 
library(ggpattern)
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

##################### Looking for specific nematode with known symbionts ####
symbionts_18s <- subset_taxa(phylo_normalized_18s, V22 == "Robbea" | V22 == "Astomonema" | V22 == "Eubostrichus")
symbionts_18s <- prune_samples(sample_sums(symbionts_18s) > 0, symbionts_18s)
symbionts_df_18s <- psmelt(symbionts_18s)
symbionts_df_18s <- subset(symbionts_df_18s, Abundance != 0)
write.csv(symbionts_df_18s, "18s_sym.csv")

symbionts_df_18s <- symbionts_df_18s %>%
  mutate(
    Type = case_when(
      SampleType %in% c("RawSediment") ~ "stripe",  # Apply pattern to these
      TRUE ~ "none"  # No pattern for others
    )
  )

symbionts_bar_18s_2 <- ggplot(
  symbionts_df_18s, 
  aes(x = Sample, y = Abundance, fill = V22, pattern = Type)
) +
  geom_bar_pattern(
    stat = "identity",
    na.rm = TRUE,
    pattern_fill = "black",       # Color of the pattern lines
    pattern_angle = 45,           # Angle of stripes
    pattern_density = 0.2,        # How dense the pattern is (0-1)
    pattern_spacing = 0.02,       # Spacing between pattern elements
    pattern_key_scale_factor = 1  # Adjust legend key size
  ) +
  facet_nested(. ~ Site + Habitat, scales = "free") +
  scale_fill_manual(values = palette, name = "Chemosynthetic Nematodes") +
  scale_pattern_manual(
    values = c("stripe" = "stripe", "none" = "none"),
               name = "Sample Type",
    labels = c("none" = "Ludox", "stripe" = "Raw Sediment")) +
  guides(fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.title.x = element_blank()) +
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold"))


palette <- c("#335c67", "#fff3b0","#e09f3e")

symbionts_bar_18s <- ggplot(symbionts_df_18s, aes(x = Sample, y = Abundance, fill = V22)) +
  geom_bar(stat = "identity", na.rm = TRUE) +
  facet_nested(. ~ Site+Habitat, scales = "free") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette, name = "Chemosynthetic Nematodes") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, vjust = 0.5, face = "bold")) + # adjusts text of x axis
  #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + 
  #theme(axis.title.x = element_blank()) +# adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold"))+  # Format facet grid title 
  scale_x_discrete(label = function(x) stringr::str_replace(x, "18S-bodega-bay_","")) 


symbionts_bar_18s

symbiont_data <- read_tsv("Raw Data/blast_hit.tsv")

custom_theme <- ttheme_minimal(
  core = list(
    bg_params = list(fill = NA, col = "black") # Ensure cell borders are black
  ),
  colhead = list(
    bg_params = list(fill = NA, col = "black") # Header cell borders
  ),
  rowhead = list(
    bg_params = list(fill = NA, col = "black") # Row header borders (if applicable)
  )
)

symbionts_plot <- grid.arrange(symbionts_bar_18s_2, tableGrob(symbiont_data, theme = custom_theme, rows = NULL), nrow = 1)

##################### Beta Diversity for Raw Sediment All Samples #################
 
# PCoA Raw Sediment
phylo_18s_raw_ord_pcoa <- ordinate(phylo_normalized_raw_18s, "PCoA", "bray")
phylo_18s_raw_pcoa <- phyloseq::plot_ordination(phylo_normalized_raw_18s, phylo_18s_raw_ord_pcoa, type="samples", color="Site", shape="Habitat", title = "18S rRNA Raw Sediment - Microeukaryotes") + 
  geom_point(size = 3) +
  scale_color_manual(values = c("#C8AB83","#2E5266","#E43F6F"),breaks = c("Campbell Cove", "Westside Park", "Mason's Marina")) +
  scale_shape_manual(values = c(1,16),
                     name = "Habitat",
                     breaks = c("Sea Grass", "Bare Sediment")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  theme(legend.position = "none")
  
phylo_18s_raw_pcoa$layers <- phylo_18s_raw_pcoa$layers[-1]

phylo_18s_raw_pcoa


# PERMANOVA for All 16s Samples
# Calculate bray curtis distance matrix
bray_18s_raw <- phyloseq::distance(phylo_normalized_raw_18s, method = "bray")

# make a data frame from the sample_data
sampledf_18s_raw <- data.frame(sample_data(phylo_normalized_raw_18s))

# Create a combined factor
sampledf_18s_raw$Site_Habitat <- interaction(sampledf_18s_raw$Site, sampledf_18s_raw$Habitat)


# Adonis test
adonis_results_18s_site_raw <- adonis2(bray_18s_raw ~ Site*Habitat, data = sampledf_18s_raw, permutations = 9999, by = "terms")
print(adonis_results_18s_site_raw) # Significance found at Pr(>F) = .001 (***)


# Homogeneity of dispersion test
beta_18s_results_site_raw <- betadisper(bray_18s_raw, sampledf_18s_raw$Site_Habitat)
plot(beta_18s_results_site_raw)
permutest(beta_18s_results_site_raw) # Significance found at Pr(>F) = .001 (***)
TukeyHSD(beta_18s_results_site_raw)


# Pairwise
results_list_18s_raw_hab <- list()
for (hab in unique(sampledf_18s_raw$Habitat)) {
  idx <- sampledf_18s_raw$Habitat == hab
  sub_dist <- as.dist(as.matrix(bray_18s_raw)[idx, idx])
  sub_code <- droplevels(as.factor(sampledf_18s_raw$Site[idx]))
  pw <- pairwise.perm.manova(sub_dist, sub_code,
                             nperm = 9999,
                             p.method = "fdr")
  results_list_18s_raw_hab[[hab]] <- pw
}
print(results_list_18s_raw_hab)


results_list_18s_raw_site <- list()
for (site in unique(sampledf_18s_raw$Site)) {
  idx <- sampledf_18s_raw$Site == site
  sub_dist <- as.dist(as.matrix(bray_18s_raw)[idx, idx])
  sub_code <- droplevels(as.factor(sampledf_18s_raw$Habitat[idx]))
  pw <- pairwise.perm.manova(sub_dist, sub_code,
                             nperm = 9999,
                             p.method = "fdr")
  results_list_18s_raw_site[[site]] <- pw
}
print(results_list_18s_raw_site)


phylo_18s_raw_pcoa_habitat <- plot_ordination(phylo_normalized_raw_18s, phylo_18s_raw_ord_pcoa, type="samples", color="Habitat", shape="Site", title = "18S rRNA Raw Sediment - Microeukaryotes") + 
  geom_point(size = 3) +
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  scale_shape_manual(values = c(15, 17, 19), breaks = c("Campbell Cove", "Westside Park", "Mason's Marina")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  theme(legend.position = "none")

phylo_18s_raw_pcoa_habitat$layers <- phylo_18s_raw_pcoa_habitat$layers[-1]

phylo_18s_raw_pcoa_habitat




# NMDS Raw Sediment
phylo_18s_raw_ord_nmds <- ordinate(phylo_normalized_raw_18s, "NMDS", "bray")
phylo_18s_raw_nmds <- plot_ordination(phylo_normalized_raw_18s, phylo_18s_raw_ord_nmds, type="samples", color="Site", shape="Habitat", title = "18S rRNA Raw Sediment - Microeukaryotes") + 
  geom_point(size = 3) +
  scale_color_manual(values = c("#C8AB83","#2E5266","#E43F6F"),breaks = c("Campbell Cove", "Westside Park", "Mason's Marina"))+
  scale_shape_manual(values = c(1,16),
                     name = "Habitat",
                     breaks = c("Sea Grass", "Bare Sediment")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = 1.1, y = 1.4, label ="2D Stress: 0.20") +
  theme(legend.position = "none")

phylo_18s_raw_nmds$layers <- phylo_18s_raw_nmds$layers[-1]

phylo_18s_raw_nmds




phylo_18s_raw_nmds_habitat <- plot_ordination(phylo_normalized_raw_18s, phylo_18s_raw_ord_nmds, type="samples", color="Habitat", shape="Site", title = "18S rRNA Raw Sediment - Microeukaryotes") + 
  geom_point(size = 3) +
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  scale_shape_manual(values = c(15, 17, 19), breaks = c("Campbell Cove", "Westside Park", "Mason's Marina")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = 1.1, y = 1.4, label ="2D Stress: 0.20") +
  theme(legend.position = "none")

phylo_18s_raw_nmds_habitat$layers <- phylo_18s_raw_nmds_habitat$layers[-1]

phylo_18s_raw_nmds_habitat 

##################### Beta Diversity for Ludox All Samples ###########################
# PCoA Ludox
#phyo_18s_ludox_no_wpbs <- subset_samples(phylo_normalized_ludox_18s, Site != "Westside Park" | Habitat != "Bare Sediment")
#psmelt(phyo_18s_ludox_no_wpbs)
phylo_18s_ludox_ord_pcoa <- ordinate(phylo_normalized_ludox_18s, "PCoA", "bray")
phylo_18s_ludox_pcoa <- phyloseq::plot_ordination(phylo_normalized_ludox_18s, phylo_18s_ludox_ord_pcoa, type="samples", color="Site", shape="Habitat", title = "18S rRNA Ludox - Meiofauna") + 
  geom_point(size = 3) +
  scale_color_manual(values = c("#C8AB83","#2E5266","#E43F6F"),breaks = c("Campbell Cove", "Westside Park", "Mason's Marina"))+
  scale_shape_manual(values = c(1,16),
                     name = "Habitat",
                     breaks = c("Sea Grass", "Bare Sediment")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  theme(legend.position = "none") +
  scale_x_reverse()

phylo_18s_ludox_pcoa$layers <- phylo_18s_ludox_pcoa$layers[-1]

phylo_18s_ludox_pcoa

# PERMANOVA for All 18s Samples
# Calculate bray curtis distance matrix
bray_18s_ludox <- phyloseq::distance(phylo_normalized_ludox_18s, method = "bray")

# make a data frame from the sample_data
sampledf_18s_ludox <- data.frame(sample_data(phylo_normalized_ludox_18s))

# Adonis test
adonis_results_18s_site_ludox <- adonis2(bray_18s_ludox ~ Site*Habitat, data = sampledf_18s_ludox, permutations = 9999, by = "terms")
adonis_results_18s_site_ludox # Significance found at Pr(>F) = .001 (***)


# Homogeneity of dispersion test
beta_18s_results_site_ludox <- betadisper(bray_18s_ludox, sampledf_18s_ludox$Site_Habitat)
plot(beta_18s_results_site_ludox)
permutest(beta_18s_results_site_ludox) # Significance found at Pr(>F) = .001 (***)
TukeyHSD(beta_18s_results_site_ludox)

# Create a combined factor
sampledf_18s_ludox$Site_Habitat <- interaction(sampledf_18s_ludox$Site, sampledf_18s_ludox$Habitat)

# Pariwise
results_list_18s_ludox_hab <- list()
for (hab in unique(sampledf_18s_ludox$Habitat)) {
  idx <- sampledf_18s_ludox$Habitat == hab
  sub_dist <- as.dist(as.matrix(bray_18s_ludox)[idx, idx])
  sub_code <- droplevels(as.factor(sampledf_18s_ludox$Site[idx]))
  pw <- pairwise.perm.manova(sub_dist, sub_code,
                             nperm = 9999,
                             p.method = "fdr")
  results_list_18s_ludox_hab[[hab]] <- pw
}
print(results_list_18s_ludox_hab)


results_list_18s_ludox_site <- list()
for (site in unique(sampledf_18s_ludox$Site)) {
  idx <- sampledf_18s_ludox$Site == site
  sub_dist <- as.dist(as.matrix(bray_18s_ludox)[idx, idx])
  sub_code <- droplevels(as.factor(sampledf_18s_ludox$Habitat[idx]))
  pw <- pairwise.perm.manova(sub_dist, sub_code,
                             nperm = 9999,
                             p.method = "fdr")
  results_list_18s_ludox_site[[site]] <- pw
}
print(results_list_18s_ludox_site)



phylo_18s_ludox_pcoa_habitat <- plot_ordination(phylo_normalized_ludox_18s, phylo_18s_ludox_ord_pcoa, type="samples", color="Habitat", shape="Site", title = "18S rRNA Ludox - Meiofauna") + 
  geom_point(size = 3) +
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  scale_shape_manual(values = c(15, 17, 19), breaks = c("Campbell Cove", "Westside Park", "Mason's Marina")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  theme(legend.position = "none")

phylo_18s_ludox_pcoa_habitat$layers <- phylo_18s_ludox_pcoa_habitat$layers[-1]

phylo_18s_ludox_pcoa_habitat




# NMDS Ludox
phylo_18s_ludox_ord_nmds <- ordinate(phylo_normalized_ludox_18s, "NMDS", "bray")
phylo_18s_ludox_nmds <- plot_ordination(phylo_normalized_ludox_18s, phylo_18s_ludox_ord_nmds, type="samples", color="Site", shape="Habitat", title = "18S rRNA Ludox - Meiofauna") + 
  geom_point(size = 3) +
  scale_color_manual(values = c("#C8AB83","#2E5266","#E43F6F"),breaks = c("Campbell Cove", "Westside Park", "Mason's Marina"))+
  scale_shape_manual(values = c(1,16),
                     name = "Habitat",
                     breaks = c("Sea Grass", "Bare Sediment")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = 1.1, y = 1, label ="2D Stress: 0.20") +
  theme(legend.position = "none")

phylo_18s_ludox_nmds$layers <- phylo_18s_ludox_nmds$layers[-1]

phylo_18s_ludox_nmds




phylo_18s_ludox_nmds_habitat <- plot_ordination(phylo_normalized_ludox_18s, phylo_18s_ludox_ord_nmds, type="samples", color="Habitat", shape="Site", title = "18S rRNA Ludox - Meiofauna") + 
  geom_point(size = 3) +
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  scale_shape_manual(values = c(15, 17, 19), breaks = c("Campbell Cove", "Westside Park", "Mason's Marina")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = 1.1, y = 1, label ="2D Stress: 0.20") +
  theme(legend.position = "none")

phylo_18s_ludox_nmds_habitat$layers <- phylo_18s_ludox_nmds_habitat$layers[-1]

phylo_18s_ludox_nmds_habitat

##################### Beta Diversity Split by Location All Samples ###########################
# Campbell Cove Only PCoA
phylo_18s_ord_pcoa_cc <- ordinate(CC_phylo_18s_normalized, "PCoA", "bray")

phylo_18s_pcoa_cc <- plot_ordination(CC_phylo_18s_normalized, phylo_18s_ord_pcoa_cc, type="samples", color="Habitat", shape = "SampleType",
                                     title = "18S rRNA Campbell Cove") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold"))

phylo_18s_pcoa_cc$layers <- phylo_18s_pcoa_cc$layers[-1]

phylo_18s_pcoa_cc


# Westside Park Only PCoA
phylo_18s_ord_pcoa_wp <- ordinate(WP_phylo_18s_normalized, "PCoA", "bray")

phylo_18s_pcoa_wp <- plot_ordination(WP_phylo_18s_normalized, phylo_18s_ord_pcoa_wp, type="samples", color="Habitat", shape = "SampleType",
                                     title = "18S rRNA Westside Park") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold"))

phylo_18s_pcoa_wp$layers <- phylo_18s_pcoa_wp$layers[-1]

phylo_18s_pcoa_wp



# Mason's Marina Only PCoA
phylo_18s_ord_pcoa_mm <- ordinate(MM_phylo_18s_normalized, "PCoA", "bray")

phylo_18s_pcoa_mm <- plot_ordination(MM_phylo_18s_normalized, phylo_18s_ord_pcoa_mm, type="samples", color="Habitat", shape = "SampleType",
                                     title = "18S rRNA Mason's Marina") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold"))

phylo_18s_pcoa_mm$layers <- phylo_18s_pcoa_mm$layers[-1]

phylo_18s_pcoa_mm




#NMDS
# Campbell Cove Only NMDS
phylo_18s_ord_nmds_cc <- ordinate(CC_phylo_18s_normalized, "NMDS", "bray")

phylo_18s_nmds_cc <- plot_ordination(CC_phylo_18s_normalized, phylo_18s_ord_nmds_cc, type="samples", color="Habitat", shape = "SampleType",
                                     title = "18S rRNA Campbell Cove") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = 1.1, y = 1.1, label ="2D Stress: 0.19")

phylo_18s_nmds_cc$layers <- phylo_18s_nmds_cc$layers[-1]

phylo_18s_nmds_cc




# Westside Park Only NMDS
phylo_18s_ord_nmds_wp <- ordinate(WP_phylo_18s_normalized, "NMDS", "bray")

phylo_18s_nmds_wp <- plot_ordination(WP_phylo_18s_normalized, phylo_18s_ord_nmds_wp, type="samples", color="Habitat", shape = "SampleType",
                                     title = "18S rRNA Westside Park") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = 1.1, y = 1.2, label ="2D Stress: 0.15")

phylo_18s_nmds_wp$layers <- phylo_18s_nmds_wp$layers[-1]

phylo_18s_nmds_wp




# Mason's Marina Only NMDS
phylo_18s_ord_nmds_mm <- ordinate(MM_phylo_18s_normalized, "NMDS", "bray")

phylo_18s_nmds_mm <- plot_ordination(MM_phylo_18s_normalized, phylo_18s_ord_nmds_mm, type="samples", color="Habitat", shape = "SampleType",
                                     title = "18S rRNA Mason's Marina") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = 1.1, y = 1.1, label ="2D Stress: 0.19")

phylo_18s_nmds_mm$layers <- phylo_18s_nmds_mm$layers[-1]

phylo_18s_nmds_mm


##################### Beta Diversity Split by Location Raw Sediment Sample Type ###########################
# Campbell Cove Only PCoA
phylo_18s_ord_pcoa_cc_raw <- ordinate(CC_phylo_18s_normalized_raw, "PCoA", "bray")

phylo_18s_pcoa_cc_raw <- plot_ordination(CC_phylo_18s_normalized_raw, phylo_18s_ord_pcoa_cc_raw, type="samples", color="Habitat", shape = "SampleType",
                                         title = "18S rRNA Raw Sediment Campbell Cove - Microeukaryotes") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold"))

phylo_18s_pcoa_cc_raw$layers <- phylo_18s_pcoa_cc_raw$layers[-1]

phylo_18s_pcoa_cc_raw

# Calculate bray curtis distance matrix
bray_18s_cc_raw <- phyloseq::distance(CC_phylo_18s_normalized_raw, method = "bray")

# make a data frame from the sample_data
sampledf_18s_cc_raw <- data.frame(sample_data(CC_phylo_18s_normalized_raw))

# Adonis test
adonis_results_18s_cc_raw <- adonis2(bray_18s_cc_raw ~ Habitat, data = sampledf_18s_cc_raw)
adonis_results_18s_cc_raw # Significance found at Pr(>F) = .001 (***)



# Westside Park Only PCoA
phylo_18s_ord_pcoa_wp_raw <- ordinate(WP_phylo_18s_normalized_raw, "PCoA", "bray")

phylo_18s_pcoa_wp_raw <- plot_ordination(WP_phylo_18s_normalized_raw, phylo_18s_ord_pcoa_wp_raw, type="samples", color="Habitat", shape = "SampleType",
                                         title = "18S rRNA Raw Sediment Westside Park - Microeukaryotes") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold"))

phylo_18s_pcoa_wp_raw$layers <- phylo_18s_pcoa_wp_raw$layers[-1]

phylo_18s_pcoa_wp_raw


# Calculate bray curtis distance matrix
bray_18s_wp_raw <- phyloseq::distance(WP_phylo_18s_normalized_raw, method = "bray")

# make a data frame from the sample_data
sampledf_18s_wp_raw <- data.frame(sample_data(WP_phylo_18s_normalized_raw))

# Adonis test
adonis_results_18s_wp_raw <- adonis2(bray_18s_wp_raw ~ Habitat, data = sampledf_18s_wp_raw)
adonis_results_18s_wp_raw # Significance found at Pr(>F) = .001 (***)


# Mason's Marina Only PCoA
phylo_18s_ord_pcoa_mm_raw <- ordinate(MM_phylo_18s_normalized_raw, "PCoA", "bray")

phylo_18s_pcoa_mm_raw <- plot_ordination(MM_phylo_18s_normalized_raw, phylo_18s_ord_pcoa_mm_raw, type="samples", color="Habitat", shape = "SampleType",
                                         title = "18S rRNA Raw Sediment Mason's Marina - Microeukaryotes") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold"))

phylo_18s_pcoa_mm_raw$layers <- phylo_18s_pcoa_mm_raw$layers[-1]

phylo_18s_pcoa_mm_raw


# Calculate bray curtis distance matrix
bray_18s_mm_raw <- phyloseq::distance(MM_phylo_18s_normalized_raw, method = "bray")

# make a data frame from the sample_data
sampledf_18s_MM_raw <- data.frame(sample_data(MM_phylo_18s_normalized_raw))

# Adonis test
adonis_results_18s_mm_raw <- adonis2(bray_18s_mm_raw ~ Habitat, data = sampledf_18s_MM_raw)
adonis_results_18s_mm_raw # Significance found at Pr(>F) = .001 (***)



#NMDS
# Campbell Cove Only NMDS
phylo_18s_ord_nmds_cc_raw <- ordinate(CC_phylo_18s_normalized_raw, "NMDS", "bray")

phylo_18s_nmds_cc_raw <- plot_ordination(CC_phylo_18s_normalized_raw, phylo_18s_ord_nmds_cc_raw, type="samples", color="Habitat", shape = "SampleType",
                                         title = "18S rRNA Raw Sediment Campbell Cove - Microeukaryotes") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = .7, y = 1.1, label ="2D Stress: 0.19")

phylo_18s_nmds_cc_raw$layers <- phylo_18s_nmds_cc_raw$layers[-1]

phylo_18s_nmds_cc_raw




# Westside Park Only NMDS
phylo_18s_ord_nmds_wp_raw <- ordinate(WP_phylo_18s_normalized_raw, "NMDS", "bray")

phylo_18s_nmds_wp_raw <- plot_ordination(WP_phylo_18s_normalized_raw, phylo_18s_ord_nmds_wp_raw, type="samples", color="Habitat", shape = "SampleType",
                                         title = "18S rRNA Raw Sediment Westside Park - Microeukaryotes") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = 1, y = 1.2, label ="2D Stress: 0.17")

phylo_18s_nmds_wp_raw$layers <- phylo_18s_nmds_wp_raw$layers[-1]

phylo_18s_nmds_wp_raw




# Mason's Marina Only NMDS
phylo_18s_ord_nmds_mm_raw <- ordinate(MM_phylo_18s_normalized_raw, "NMDS", "bray")

phylo_18s_nmds_mm_raw <- plot_ordination(MM_phylo_18s_normalized_raw, phylo_18s_ord_nmds_mm_raw, type="samples", color="Habitat", shape = "SampleType",
                                         title = "18S rRNA Raw Sediment Mason's Marina - Microeukaryotes") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = .8, y = .7, label ="2D Stress: 0.15")

phylo_18s_nmds_mm_raw$layers <- phylo_18s_nmds_mm_raw$layers[-1]

phylo_18s_nmds_mm_raw




##################### Beta Diversity Split by Location Ludox Sample Type ###########################
# Campbell Cove Only PCoA
phylo_18s_ord_pcoa_cc_ludox <- ordinate(CC_phylo_18s_normalized_ludox, "PCoA", "bray")

phylo_18s_pcoa_cc_ludox <- plot_ordination(CC_phylo_18s_normalized_ludox, phylo_18s_ord_pcoa_cc_ludox, type="samples", color="Habitat", shape = "SampleType",
                                           title = "18S rRNA Ludox Campbell Cove - Meiofauna") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold"))

phylo_18s_pcoa_cc_ludox$layers <- phylo_18s_pcoa_cc_ludox$layers[-1]

phylo_18s_pcoa_cc_ludox

# Calculate bray curtis distance matrix
bray_18s_cc_ludox <- phyloseq::distance(CC_phylo_18s_normalized_ludox, method = "bray")

# make a data frame from the CC_phylo_18s_normalized_ludox
sampledf_18s_cc_ludox <- data.frame(sample_data(CC_phylo_18s_normalized_ludox))

# Adonis test
adonis_results_18s_cc_ludox <- adonis2(bray_18s_cc_ludox ~ Habitat, data = sampledf_18s_cc_ludox)
adonis_results_18s_cc_ludox # Significance found at Pr(>F) = .001 (***)


# Westside Park Only PCoA
phylo_18s_ord_pcoa_wp_ludox <- ordinate(WP_phylo_18s_normalized_ludox, "PCoA", "bray")

phylo_18s_pcoa_wp_ludox <- plot_ordination(WP_phylo_18s_normalized_ludox, phylo_18s_ord_pcoa_wp_ludox, type="samples", color="Habitat", shape = "SampleType",
                                           title = "18S rRNA Ludox Westside Park - Meiofauna") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold"))

phylo_18s_pcoa_wp_ludox$layers <- phylo_18s_pcoa_wp_ludox$layers[-1]

phylo_18s_pcoa_wp_ludox

# Calculate bray curtis distance matrix
bray_18s_wp_ludox <- phyloseq::distance(WP_phylo_18s_normalized_ludox, method = "bray")

# make a data frame from the CC_phylo_18s_normalized_ludox
sampledf_18s_wp_ludox <- data.frame(sample_data(WP_phylo_18s_normalized_ludox))

# Adonis test
adonis_results_18s_wp_ludox <- adonis2(bray_18s_wp_ludox ~ Habitat, data = sampledf_18s_wp_ludox)
adonis_results_18s_wp_ludox # Significance found at Pr(>F) = .001 (***)

# Mason's Marina Only PCoA
phylo_18s_ord_pcoa_mm_ludox <- ordinate(MM_phylo_18s_normalized_ludox, "PCoA", "bray")

phylo_18s_pcoa_mm_ludox <- plot_ordination(MM_phylo_18s_normalized_ludox, phylo_18s_ord_pcoa_mm_ludox, type="samples", color="Habitat", shape = "SampleType",
                                           title = "18S rRNA Ludox Mason's Marina - Meiofauna") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold"))

phylo_18s_pcoa_mm_ludox$layers <- phylo_18s_pcoa_mm_ludox$layers[-1]

phylo_18s_pcoa_mm_ludox


# Calculate bray curtis distance matrix
bray_18s_mm_ludox <- phyloseq::distance(MM_phylo_18s_normalized_ludox, method = "bray")

# make a data frame from the CC_phylo_18s_normalized_ludox
sampledf_18s_mm_ludox <- data.frame(sample_data(MM_phylo_18s_normalized_ludox))

# Adonis test
adonis_results_18s_mm_ludox <- adonis2(bray_18s_mm_ludox ~ Habitat, data = sampledf_18s_mm_ludox)
adonis_results_18s_mm_ludox # Significance found at Pr(>F) = .002 (**)



#NMDS
# Campbell Cove Only NMDS
phylo_18s_ord_nmds_cc_ludox <- ordinate(CC_phylo_18s_normalized_ludox, "NMDS", "bray")

phylo_18s_nmds_cc_ludox <- plot_ordination(CC_phylo_18s_normalized_ludox, phylo_18s_ord_nmds_cc_ludox, type="samples", color="Habitat", shape = "SampleType",
                                           title = "18S rRNA Ludox Campbell Cove - Meiofauna") +
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = 1, y = .7, label ="2D Stress: 0.20")

phylo_18s_nmds_cc_ludox$layers <- phylo_18s_nmds_cc_ludox$layers[-1]

phylo_18s_nmds_cc_ludox




# Westside Park Only NMDS
phylo_18s_ord_nmds_wp_ludox <- ordinate(WP_phylo_18s_normalized_ludox, "NMDS", "bray")

phylo_18s_nmds_wp_ludox <- plot_ordination(WP_phylo_18s_normalized_ludox, phylo_18s_ord_nmds_wp_ludox, type="samples", color="Habitat", shape = "SampleType",
                                           title = "18S rRNA Ludox Westside Park - Meiofauna") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = 0, y = 1, label ="2D Stress: 0.11")

phylo_18s_nmds_wp_ludox$layers <- phylo_18s_nmds_wp_ludox$layers[-1]

phylo_18s_nmds_wp_ludox




# Mason's Marina Only NMDS
phylo_18s_ord_nmds_mm_ludox <- ordinate(MM_phylo_18s_normalized_ludox, "NMDS", "bray")

phylo_18s_nmds_mm_ludox <- plot_ordination(MM_phylo_18s_normalized_ludox, phylo_18s_ord_nmds_mm_ludox, type="samples", color="Habitat", shape = "SampleType",
                                           title = "18S rRNA Ludox Mason's Marina - Meiofauna") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = 1, y = 1, label ="2D Stress: 0.20")

phylo_18s_nmds_mm_ludox$layers <- phylo_18s_nmds_mm_ludox$layers[-1]

phylo_18s_nmds_mm_ludox




##################### Beta Diversity Nematodes Only ####################
# PCoA 
phylo_18s_nematode_ord_pcoa <- ordinate(phylo_normalized_nematoda_18s, "PCoA", "bray")
phylo_18s_nematode_pcoa <- plot_ordination(phylo_normalized_nematoda_18s, phylo_18s_nematode_ord_pcoa, type="samples", color="Site", shape="Habitat",
                                           title = "18S rRNA - Nematodes") + 
  geom_point(size = 3) +
  scale_color_manual(values = c("#C8AB83","#2E5266","#E43F6F"),breaks = c("Campbell Cove", "Westside Park", "Mason's Marina"))+
  scale_shape_manual(values = c(1,16),
                     name = "Habitat",
                     breaks = c("Sea Grass", "Bare Sediment")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold"))

phylo_18s_nematode_pcoa$layers <- phylo_18s_nematode_pcoa$layers[-1]

phylo_18s_nematode_pcoa




phylo_18s_nematode_pcoa_habitat <- plot_ordination(phylo_normalized_nematoda_18s, phylo_18s_nematode_ord_pcoa, type="samples", color="Habitat", shape="Site",
                                                   title = "18S rRNA - Nematodes") + 
  geom_point(size = 3) +
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  scale_shape_manual(values = c(15, 17, 19), breaks = c("Campbell Cove", "Westside Park", "Mason's Marina")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold"))

phylo_18s_nematode_pcoa_habitat$layers <- phylo_18s_nematode_pcoa_habitat$layers[-1]

phylo_18s_nematode_pcoa_habitat





# NMDS Nematode
phylo_18s_nematode_ord_nmds <- ordinate(phylo_normalized_nematoda_18s, "NMDS", "bray")
phylo_18s_nematode_nmds <- plot_ordination(phylo_normalized_nematoda_18s, phylo_18s_nematode_ord_nmds, type="samples", color="Site", shape="Habitat",
                                           title = "18S rRNA - Nematodes") + 
  geom_point(size = 3) +
  scale_color_manual(values = c("#C8AB83","#2E5266","#E43F6F"),breaks = c("Campbell Cove", "Westside Park", "Mason's Marina"))+
  scale_shape_manual(values = c(1,16),
                     name = "Habitat",
                     breaks = c("Sea Grass", "Bare Sediment")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = 1, y = 1, label ="2D Stress: 0.20")

phylo_18s_nematode_nmds$layers <- phylo_18s_nematode_nmds$layers[-1]

phylo_18s_nematode_nmds




phylo_18s_nematode_nmds_habitat <- plot_ordination(phylo_normalized_nematoda_18s, phylo_18s_nematode_ord_nmds, type="samples", color="Habitat", shape="Site",
                                                   title = "18S rRNA - Nematodes") + 
  geom_point(size = 3) +
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  scale_shape_manual(values = c(15, 17, 19), breaks = c("Campbell Cove", "Westside Park", "Mason's Marina")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = 1, y = 1, label ="2D Stress: 0.20")

phylo_18s_nematode_nmds_habitat$layers <- phylo_18s_nematode_nmds_habitat$layers[-1]

phylo_18s_nematode_nmds_habitat





# Now by Location
# Campbell Cove Only PCoA
phylo_18s_ord_pcoa_cc_nematode <- ordinate(CC_phylo_18s_normalized_nematode, "PCoA", "bray")

phylo_18s_pcoa_cc_nematode <- plot_ordination(CC_phylo_18s_normalized_nematode, phylo_18s_ord_pcoa_cc_nematode, type="samples", color="Habitat",
                                              title = "18S rRNA Campbell Cove - Nematodes") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold"))

phylo_18s_pcoa_cc_nematode$layers <- phylo_18s_pcoa_cc_nematode$layers[-1]

phylo_18s_pcoa_cc_nematode


# Westside Park Only PCoA
phylo_18s_ord_pcoa_wp_nematode <- ordinate(WP_phylo_18s_normalized_nematode, "PCoA", "bray")

phylo_18s_pcoa_wp_nematode <- plot_ordination(WP_phylo_18s_normalized_nematode, phylo_18s_ord_pcoa_wp_nematode, type="samples", color="Habitat",
                                              title = "18S rRNA Westside Park - Nematodes") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold"))

phylo_18s_pcoa_wp_nematode$layers <- phylo_18s_pcoa_wp_nematode$layers[-1]

phylo_18s_pcoa_wp_nematode



# Mason's Marina Only PCoA
phylo_18s_ord_pcoa_mm_nematode <- ordinate(MM_phylo_18s_normalized_nematode, "PCoA", "bray")

phylo_18s_pcoa_mm_nematode <- plot_ordination(MM_phylo_18s_normalized_nematode, phylo_18s_ord_pcoa_mm_nematode, type="samples", color="Habitat",
                                              title = "18S rRNA Mason's Marina - Nematodes") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold"))

phylo_18s_pcoa_mm_nematode$layers <- phylo_18s_pcoa_mm_nematode$layers[-1]

phylo_18s_pcoa_mm_nematode





#NMDS
# Campbell Cove Only NMDS
phylo_18s_ord_nmds_cc_nematode <- ordinate(CC_phylo_18s_normalized_nematode, "NMDS", "bray")

phylo_18s_nmds_cc_nematode <- plot_ordination(CC_phylo_18s_normalized_nematode, phylo_18s_ord_nmds_cc_nematode, type="samples", color="Habitat",
                                              title = "18S rRNA Campbell Cove - Nematodes") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = 1.5, y = 1, label ="2D Stress: 0.17")

phylo_18s_nmds_cc_nematode$layers <- phylo_18s_nmds_cc_nematode$layers[-1]

phylo_18s_nmds_cc_nematode




# Westside Park Only NMDS
phylo_18s_ord_nmds_wp_nematode <- ordinate(WP_phylo_18s_normalized_nematode, "NMDS", "bray")

phylo_18s_nmds_wp_nematode <- plot_ordination(WP_phylo_18s_normalized_nematode, phylo_18s_ord_nmds_wp_nematode, type="samples", color="Habitat",
                                              title = "18S rRNA Westside Park - Nematodes") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = -1, y = 1, label ="2D Stress: 0.12")

phylo_18s_nmds_wp_nematode$layers <- phylo_18s_nmds_wp_nematode$layers[-1]

phylo_18s_nmds_wp_nematode




# Mason's Marina Only NMDS
phylo_18s_ord_nmds_mm_nematode <- ordinate(MM_phylo_18s_normalized_nematode, "NMDS", "bray")

phylo_18s_nmds_mm_nematode <- plot_ordination(MM_phylo_18s_normalized_nematode, phylo_18s_ord_nmds_mm_nematode, type="samples", color="Habitat",
                                              title = "18S rRNA Mason's Marina - Nematodes") + 
  scale_color_manual(values = c("saddlebrown","#00A572"),breaks = c("Bare Sediment", "Sea Grass")) +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  geom_point(size = 3) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(color = "black", face = "bold")) +
  annotate("text", x = .8, y = .8, label ="2D Stress: 0.22")

phylo_18s_nmds_mm_nematode$layers <- phylo_18s_nmds_mm_nematode$layers[-1]

phylo_18s_nmds_mm_nematode

##################### Nematode taxonomy bar charts Top 10 and Top 20 at Family and Genus Levels ###################
# Family
# Function to collapse a certain number of taxa into category others
merge_top10_18s_nema <- function(phylo_nematoda_18s, top=9){
  transformed <- transform_sample_counts(phylo_nematoda_18s, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:22 if there are species level
  }
  return(merged)
}

# Agglomerated taxa down to V21 rank 
# Use tax_fix function to assign taxa names to lower ranks that are unknown
temp_phylo_18s_nema <- phylo_nematoda_18s

fix_phylo_18s_nema <- tax_fix(temp_phylo_18s_nema, unknowns = c("Unassigned", "uncultured", "Unkown", "uncultured eukaryote", "Unknown Family"))

glom_18s_nema <- tax_glom(fix_phylo_18s_nema, taxrank = "V21")

# Run function on phyloseq object and make into data frame
phy_18_nema_top10_V21 <- merge_top10_18s_nema(glom_18s_nema, top=9)

get_taxa_unique(phy_18_nema_top10_V21, taxonomic.rank = "V21")

phy_18_nema_top10_V21_df <- psmelt(phy_18_nema_top10_V21)
unique(phy_18_nema_top10_V21_df$V21)

phy_18_nema_top10_V21_df$V21 <- phy_18_nema_top10_V21_df$V21 %>%
  replace_na('Others')


# Add common factors to use for plotting
phy_18_nema_top10_agr = aggregate(Abundance~Sample+Site+Habitat+V21, data=phy_18_nema_top10_V21_df, FUN=mean) 
unique(phy_18_nema_top10_agr$V21)

# Put "Others" to the final of the Family list - top 10
phy_18_nema_top10_agr$V21 <- factor(phy_18_nema_top10_agr$V21,
                                   levels = c("Anticomidae", "Chromadoridae", "Comesomatidae",
                                              "Cyatholaimidae", "Enchelidiidae", "Enoplida V17",
                                              "Linhomoeidae", "Nematoda V14", "Oncholaimidae",
                                              "Xyalidae", "Others"))

# Reorder Site levels
phy_18_nema_top10_agr$Site = factor(phy_18_nema_top10_agr$Site, levels=c("Campbell Cove","Westside Park","Mason's Marina"))

# Color palette
colors_top10_3 <- c("#67001F", "#B2182B" ,"#D6604D" ,"#F4A582" ,"#FDDBC7" ,"#D1E5F0", "#92C5DE", "#4393C3", "#2166AC","grey")

# Plot by site - Family level top 10
taxonomy_bar_18s_nema_fam <- ggplot(phy_18_nema_top10_agr, aes(x = Sample, y = Abundance, fill = V21)) +
  facet_nested(. ~ Site+Habitat, scales = "free",
               labeller = labeller(Site = site.labs, Habitat = habitat.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  #geom_text(aes(label = ifelse(round(Abundance*100) >= 1, paste(round(Abundance*100, digits = 0), "%"), "")), size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top10_4, name = "Family") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
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

taxonomy_bar_18s_nema_fam

class(phy_18_nema_top10_agr)
write.csv(phy_18_nema_top10_agr, "phy_18_nema_top10_agr.csv")


# Genus
# Agglomerated taxa down to V22 rank 
glom_18s_nema_genus <- tax_glom(fix_phylo_18s_nema, taxrank = "V22")

# Run function on phyloseq object and make into data frame
phy_18_nema_top10_V22 <- merge_top10_18s_nema(glom_18s_nema_genus, top=9)

get_taxa_unique(phy_18_nema_top10_V22, taxonomic.rank = "V22")

phy_18_nema_top10_V22_df <- psmelt(phy_18_nema_top10_V22)
unique(phy_18_nema_top10_V22_df$V22)

phy_18_nema_top10_V22_df$V22 <- phy_18_nema_top10_V22_df$V22 %>%
  replace_na('Others')


# Add common factors to use for plotting
phy_18_nema_top10_agr_V22 = aggregate(Abundance~Sample+Site+Habitat+V22, data=phy_18_nema_top10_V22_df, FUN=mean) 
unique(phy_18_nema_top10_agr_V22$V22)

# Put "Others" to the final of the Family list - top 10
phy_18_nema_top10_agr_V22$V22 <- factor(phy_18_nema_top10_agr_V22$V22,
                                    levels = c("Anticoma", "Calyptronema", "Daptonema",
                                               "Enoplida V17", "Metalinhomoeus", "Ptycholaimellus",
                                               "Sabatieria", "Terschellingia", "Viscosia", "Others"))

# Reorder Site levels
phy_18_nema_top10_agr_V22$Site = factor(phy_18_nema_top10_agr_V22$Site, levels=c("Campbell Cove","Westside Park","Mason's Marina"))

# Plot by site - Genus level top 10
taxonomy_bar_18s_nema_genus <- ggplot(phy_18_nema_top10_agr_V22, aes(x = Sample, y = Abundance, fill = V22)) +
  facet_nested(. ~ Site+Habitat, scales = "free",
               labeller = labeller(Site = site.labs, Habitat = habitat.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  #geom_text(aes(label = ifelse(round(Abundance*100) >= 1, paste(round(Abundance*100, digits = 0), "%"), "")), size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top10_4, name = "Genus") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, vjust = 0.5, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
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

taxonomy_bar_18s_nema_genus


write.csv(phy_18_nema_top10_agr_V22, "phy_18_nema_top10_agr_V22.csv")





# Top 20
# Family
# Function to collapse a certain number of taxa into category others
merge_top20_18s_nema <- function(phylo_dataframe, top=19){
  transformed <- transform_sample_counts(phylo_dataframe, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:23] <- "Others"} # 1:22 if there are species level
  }
  return(merged)
}


# Run function on phyloseq object and make into data frame
phy_18_nema_top20_V21 <- merge_top20_18s_nema(glom_18s_nema, top=19)

get_taxa_unique(phy_18_nema_top20_V21, taxonomic.rank = "V21")

phy_18_nema_top20_V21_df <- psmelt(phy_18_nema_top20_V21)
unique(phy_18_nema_top20_V21_df$V21)

phy_18_nema_top20_V21_df$V21 <- phy_18_nema_top20_V21_df$V21 %>%
  replace_na('Others')


# Add common factors to use for plotting
phy_18_nema_top20_agr = aggregate(Abundance~Sample+Site+Habitat+V21, data=phy_18_nema_top20_V21_df, FUN=mean) 
unique(phy_18_nema_top20_agr$V21)

# Put "Others" to the final of the Family list - top 20
phy_18_nema_top20_agr$V21 <- factor(phy_18_nema_top20_agr$V21,
                                    levels = c("Anticomidae", "Camacolaimidae", "Chromadoridae","Comesomatidae","Cyatholaimidae", "Desmodoridae", "Diplopeltidae",
                                               "Enchelidiidae","Enoplida V17", "Leptolaimidae", "Linhomoeidae", "Microlaimidae", "Monhysterida V17", "Monhysteridae",
                                               "Oncholaimidae","Oxystominidae","Siphonolaimidae","Tripyloididae", "Xyalidae", "Others"))

colors_top20_3 <- c("#dd3497", "#ae017e","#7a0177","#3690c0", "#74a9cf", "#000075", "#a6bddb", "#d0d1e6",
                    "#014636", "#016c59", "#02818a", "#41b6c4", "#7fcdbb","#c7e9b4","#e0f3db", "#ccece6",  "#f768a1", "#fa9fb5", "#fcc5c0","gray")

# Reorder Site levels
phy_18_nema_top20_agr$Site = factor(phy_18_nema_top10_agr$Site, levels=c("Campbell Cove","Westside Park","Mason's Marina"))

# Plot by site - Family level top 20
taxonomy_bar_18s_nema_fam_top20 <- ggplot(phy_18_nema_top20_agr, aes(x = Sample, y = Abundance, fill = V21)) +
  facet_nested(. ~ Site+Habitat, scales = "free",
               labeller = labeller(Site = site.labs, Habitat = habitat.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  #geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20_3, name = "18S Nematode Family") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_blank()) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
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

taxonomy_bar_18s_nema_fam_top20





# Run function on phyloseq object and make into data frame
phy_18_nema_top20_V22 <- merge_top20_18s_nema(glom_18s_nema_genus, top=19)

get_taxa_unique(phy_18_nema_top20_V22, taxonomic.rank = "V22")

phy_18_nema_top20_V22_df <- psmelt(phy_18_nema_top20_V22)
unique(phy_18_nema_top20_V22_df$V22)

phy_18_nema_top20_V22_df$V22 <- phy_18_nema_top20_V22_df$V22 %>%
  replace_na('Others')


# Add common factors to use for plotting
phy_18_nema_top20_agr_V22 = aggregate(Abundance~Sample+Site+Habitat+V22, data=phy_18_nema_top20_V22_df, FUN=mean) 
unique(phy_18_nema_top20_agr_V22$V22)

# Put "Others" to the final of the Family list - top 10
phy_18_nema_top20_agr_V22$V22 <- factor(phy_18_nema_top20_agr_V22$V22,
                                        levels = c("Acanthonchus", "Anticoma", "Calyptronema", "Chromadoropsis", "Daptonema", "Deontolaimus", "Desmolaimus",
                                                   "Dichromadora", "Enoplida V17", "Leptolaimus", "Metalinhomoeus", "Microlaimus", "Monhysterida V17",
                                                   "Paracanthonchus", "Ptycholaimellus", "Sabatieria", "Spilophorella", "Terschellingia", "Viscosia", "Others"))
                                          
# Reorder Site levels
phy_18_nema_top20_agr_V22$Site = factor(phy_18_nema_top20_agr_V22$Site, levels=c("Campbell Cove","Westside Park","Mason's Marina"))

# Plot by site - Genus level top 20
taxonomy_bar_18s_nema_genus_top20 <- ggplot(phy_18_nema_top20_agr_V22, aes(x = Sample, y = Abundance, fill = V22)) +
  facet_nested(. ~ Site+Habitat, scales = "free",
               labeller = labeller(Site = site.labs, Habitat = habitat.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  #geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20_3, name = "18S Nematode Genus") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, vjust = 0.5, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12))+ # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
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

taxonomy_bar_18s_nema_genus_top20

##################### Alpha Div #######
# Calculate alpha-diversity measures (For plot purposes only!)
# This can be done using the different phyloseq alpha diversity measures
# You will get a Warning message for each index since there is no singletons on the dataset

alpha_div_18S_ludox <- data.frame(
  "Observed" = phyloseq::estimate_richness(phylo_ludox_18s, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phylo_ludox_18s, measures = "Shannon"),
  "InvSimpson" = phyloseq::estimate_richness(phylo_ludox_18s, measures = "InvSimpson"),
  "Site" = phyloseq::sample_data(phylo_ludox_18s)$Site,
  "Habitat" = phyloseq::sample_data(phylo_ludox_18s)$Habitat)

alpha_div_18S_ludox$Evenness <- alpha_div_18S_ludox$Shannon/log(alpha_div_18S_ludox$Observed)

head(alpha_div_18S_ludox)

# Rename variable InvSimpson to Simpson
# The function rename & %>% works on dplyr. make sure it is loaded.
alpha_div_18S_ludox <- alpha_div_18S_ludox %>%
  dplyr::rename(Simpson = InvSimpson)
head(alpha_div_18S_ludox)

# Reorder dataframe, first categorical then numerical variables  
alpha_div_18S_ludox <- alpha_div_18S_ludox[, c( 5, 4, 1, 2, 3, 6)]
head(alpha_div_18S_ludox)

#Summarize alpha diversity measures by location
summary_alpha_18S_ludox_site <- alpha_div_18S_ludox %>%
  group_by(Site) %>%
  dplyr::summarise(count = n(),
                   mean_observed = mean(Observed),
                   sd_observed = sd(Observed),
                   mean_shannon = mean(Shannon),
                   sd_shannon = sd(Shannon),
                   mean_Simpson = mean(Simpson),
                   sd_Simpson = sd(Simpson),
                   mean_evenness = mean(Evenness),
                   sd_evenness = sd(Evenness))

write_csv(summary_alpha_18S_ludox_site, "summary_alpha_18S_ludox_site.csv")

#Summarize alpha diversity measures by habitat
summary_alpha_18S_ludox_habitat <- alpha_div_18S_ludox %>%
  group_by(Habitat) %>%
  dplyr::summarise(count = n(),
                   mean_observed = mean(Observed),
                   sd_observed = sd(Observed),
                   mean_shannon = mean(Shannon),
                   sd_shannon = sd(Shannon),
                   mean_Simpson = mean(Simpson),
                   sd_Simpson = sd(Simpson),
                   mean_evenness = mean(Evenness),
                   sd_evenness = sd(Evenness))

write_csv(summary_alpha_18S_ludox_habitat, "summary_alpha_18S_ludox_habitat.csv")

#Summarize alpha diversity measures by habitat within site
summary_alpha_18S_ludox_site_habitat <- alpha_div_18S_ludox %>%
  group_by(Site, Habitat) %>%
  dplyr::summarise(count = n(),
                   mean_observed = mean(Observed),
                   median_observed = median(Observed),
                   sd_observed = sd(Observed),
                   mean_shannon = mean(Shannon),
                   sd_shannon = sd(Shannon),
                   mean_Simpson = mean(Simpson),
                   sd_Simpson = sd(Simpson),
                   mean_evenness = mean(Evenness),
                   sd_evenness = sd(Evenness))

write_csv(summary_alpha_18S_ludox_site_habitat, "summary_alpha_18S_ludox_site_habitat.csv")


# KW analysis alpha_16S_nem on all metrics at once
# Remember, numerical variables are from columns 3-5. The test is by location, column 2
kw_alpha_18S_ludox <- as.data.frame(sapply(3:6, function(x) kruskal.test(alpha_div_18S_ludox[,x],
                                                                       alpha_div_18S_ludox[,2])))

# Rename columns with the proper variable names
kw_alpha_18S_ludox <- kw_alpha_18S_ludox %>%
  dplyr::rename(Observed = V1,
         Shannon = V2,
         Simpson = V3,
         Evenness = V4)

kw_alpha_18S_ludox <- t(kw_alpha_18S_ludox) # transpose
kw_alpha_18S_ludox <- as_tibble(kw_alpha_18S_ludox, rownames = "Metric") # adding rownames as a column
kw_alpha_18S_ludox <- kw_alpha_18S_ludox[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_18S_ludox) # checking object class

# KW analysis alpha_16S_nem on all metrics at once
# Remember, numerical variables are from columns 3-5. The test is by habitat, column 2
kw_alpha_18S_ludox_habitat <- as.data.frame(sapply(3:6, function(x) kruskal.test(alpha_div_18S_ludox[,x],
                                                                         alpha_div_18S_ludox[,1])))

# Rename columns with the proper variable names
kw_alpha_18S_ludox_habitat <- kw_alpha_18S_ludox_habitat %>%
  dplyr::rename(Observed = V1,
                Shannon = V2,
                Simpson = V3,
                Evenness = V4)

kw_alpha_18S_ludox_habitat <- t(kw_alpha_18S_ludox_habitat) # transpose
kw_alpha_18S_ludox_habitat <- as_tibble(kw_alpha_18S_ludox_habitat, rownames = "Metric") # adding rownames as a column
kw_alpha_18S_ludox_habitat <- kw_alpha_18S_ludox_habitat[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_18S_ludox_habitat) # checking object class

# Remember, numerical variables are from columns 3-5. The test is by site x habitat, column 7
alpha_div_18S_ludox_with_new_factor <- alpha_div_18S_ludox

alpha_div_18S_ludox_with_new_factor$Site_Habitat <- interaction(alpha_div_18S_ludox_with_new_factor$Site, alpha_div_18S_ludox_with_new_factor$Habitat)

kw_alpha_18S_ludox_sitexhabitat <- as.data.frame(sapply(3:6, function(x) kruskal.test(alpha_div_18S_ludox_with_new_factor[,x],
                                                                                 alpha_div_18S_ludox_with_new_factor[,7])))

# Rename columns with the proper variable names
kw_alpha_18S_ludox_sitexhabitat <- kw_alpha_18S_ludox_sitexhabitat %>%
  dplyr::rename(Observed = V1,
                Shannon = V2,
                Simpson = V3,
                Evenness = V4)

kw_alpha_18S_ludox_sitexhabitat <- t(kw_alpha_18S_ludox_sitexhabitat) # transpose
kw_alpha_18S_ludox_sitexhabitat <- as_tibble(kw_alpha_18S_ludox_sitexhabitat, rownames = "Metric") # adding rownames as a column
kw_alpha_18S_ludox_sitexhabitat <- kw_alpha_18S_ludox_sitexhabitat[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_18S_ludox_sitexhabitat) # checking object class

# Save resulting table with fwrite to avoid any issues with characters
data.table::fwrite(kw_alpha_18S_ludox, "kw_alpha_18S_ludox.csv")
data.table::fwrite(kw_alpha_18S_ludox_habitat, "kw_alpha_18S_ludox_habitat.csv")
data.table::fwrite(kw_alpha_18S_ludox_sitexhabitat, "kw_alpha_18S_ludox_sitexhabitat.csv")

# Plot alpha diversity measures
# Change names
alpha_div_18S_ludox$Site = factor(alpha_div_18S_ludox$Site, levels=c("Campbell Cove","Westside Park","Mason's Marina"))

site_comparisons <- list( c("Campbell Cove", "Westside Park"), c("Westside Park", "Mason's Marina"), c("Campbell Cove", "Mason's Marina"))

alpha_color_site <- c("#C8AB83","#2E5266","#E43F6F")

ad_18s_ludox_site <- alpha_div_18S_ludox %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Simpson", "Evenness")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Simpson", "Evenness"))) %>%
  ggplot(aes(x = Site, y = value)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  stat_compare_means() +
  stat_compare_means(comparisons = site_comparisons, p.adjust.methods = "BH", aes(label = ..p.signif..), size = 4, hide.ns = FALSE) +
  geom_jitter(aes(color = Site), height = 0, width = .2) +
  facet_nested(metric ~ Habitat, scales = "free") +
  scale_color_manual(values = alpha_color_site) +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_blank()) + # adjusts the title of x axis +
  scale_y_continuous(expand = c(0.1, 0.1)) +
  scale_x_discrete(
    labels = c("CC", "WP", "MM"),
    expand = c(0.2, 0.2),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) +  # adjusts the title of the legend
  ylab("") + # add the title on y axis
  xlab("Site") +  # add the title on x axis
  theme(legend.position="none")  # add the title on x axis

ad_18s_ludox_site


alpha_color_habitat <- c("saddlebrown","#00A572")

habitat_comparisons <- list(c("Bare Sediment", "Sea Grass"))

ad_18s_ludox_habitat_ <- alpha_div_18S_ludox %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Simpson", "Evenness")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Simpson", "Evenness"))) %>%
  ggplot(aes(x = Habitat, y = value)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  stat_compare_means() +
  stat_compare_means(comparisons = habitat_comparisons, p.adjust.methods = "BH", aes(label =..p.signif..), size = 4, hide.ns = FALSE) +
  geom_jitter(aes(color = Habitat), height = 0, width = .2) +
  facet_nested(metric ~ Site, scales = "free") +
  scale_color_manual(values = alpha_color_habitat) +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_blank()) + # adjusts the title of x axis 
  scale_y_continuous(expand = c(0.1, 0.1)) +
  scale_x_discrete(
    labels = c("BS", "SG"),
    expand = c(0.2, 0.2),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) +  # adjusts the title of the legend
  ylab("") + # add the title on y axis
  xlab("Site") +  # add the title on x axis
  theme(legend.position="none") 

ad_18s_ludox_habitat






#### Raw
alpha_div_18S_raw <- data.frame(
  "Observed" = phyloseq::estimate_richness(phylo_raw_18s, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phylo_raw_18s, measures = "Shannon"),
  "InvSimpson" = phyloseq::estimate_richness(phylo_raw_18s, measures = "InvSimpson"),
  "Site" = phyloseq::sample_data(phylo_raw_18s)$Site,
  "Habitat" = phyloseq::sample_data(phylo_raw_18s)$Habitat)

alpha_div_18S_raw$Evenness <- alpha_div_18S_raw$Shannon/log(alpha_div_18S_raw$Observed)

head(alpha_div_18S_raw)

# Rename variable InvSimpson to Simpson
# The function rename & %>% works on dplyr. make sure it is loaded.
alpha_div_18S_raw <- alpha_div_18S_raw %>%
  dplyr::rename(Simpson = InvSimpson)
head(alpha_div_18S_raw)

# Reorder dataframe, first categorical then numerical variables  
alpha_div_18S_raw <- alpha_div_18S_raw[, c( 5, 4, 1, 2, 3, 6)]
head(alpha_div_18S_raw)

#Summarize alpha diversity measures by location
summary_alpha_18S_raw_site <- alpha_div_18S_raw %>%
  group_by(Site) %>%
  dplyr::summarise(count = n(),
                   mean_observed = mean(Observed),
                   sd_observed = sd(Observed),
                   mean_shannon = mean(Shannon),
                   sd_shannon = sd(Shannon),
                   mean_Simpson = mean(Simpson),
                   sd_Simpson = sd(Simpson),
                   mean_evenness = mean(Evenness),
                   sd_evenness = sd(Evenness))

write_csv(summary_alpha_18S_raw_site, "summary_alpha_18S_raw_site.csv")

#Summarize alpha diversity measures by habitat
summary_alpha_18S_raw_habitat <- alpha_div_18S_raw %>%
  group_by(Habitat) %>%
  dplyr::summarise(count = n(),
                   mean_observed = mean(Observed),
                   sd_observed = sd(Observed),
                   mean_shannon = mean(Shannon),
                   sd_shannon = sd(Shannon),
                   mean_Simpson = mean(Simpson),
                   sd_Simpson = sd(Simpson),
                   mean_evenness = mean(Evenness),
                   sd_evenness = sd(Evenness))

write_csv(summary_alpha_18S_raw_habitat, "summary_alpha_18S_raw_habitat.csv")

#Summarize alpha diversity measures by habitat within site
summary_alpha_18S_raw_site_habitat <- alpha_div_18S_raw %>%
  group_by(Site, Habitat) %>%
  dplyr::summarise(count = n(),
                   mean_observed = mean(Observed),
                   median_observed = median(Observed),
                   sd_observed = sd(Observed),
                   mean_shannon = mean(Shannon),
                   sd_shannon = sd(Shannon),
                   mean_Simpson = mean(Simpson),
                   sd_Simpson = sd(Simpson),
                   mean_evenness = mean(Evenness),
                   sd_evenness = sd(Evenness))

write_csv(summary_alpha_18S_raw_site_habitat, "summary_alpha_18S_raw_site_habitat.csv")


# KW analysis alpha_16S_nem on all metrics at once
# Remember, numerical variables are from columns 3-5. The test is by location, column 2
kw_alpha_18S_raw <- as.data.frame(sapply(3:6, function(x) kruskal.test(alpha_div_18S_raw[,x],
                                                                       alpha_div_18S_raw[,2])))

# Rename columns with the proper variable names
kw_alpha_18S_raw <- kw_alpha_18S_raw %>%
  dplyr::rename(Observed = V1,
         Shannon = V2,
         Simpson = V3,
         Evenness = V4)

kw_alpha_18S_raw <- t(kw_alpha_18S_raw) # transpose
kw_alpha_18S_raw <- as_tibble(kw_alpha_18S_raw, rownames = "Metric") # adding rownames as a column
kw_alpha_18S_raw <- kw_alpha_18S_raw[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_18S_raw) # checking object class


# KW analysis alpha_16S_nem on all metrics at once
# Remember, numerical variables are from columns 3-5. The test is by hab, column 1
kw_alpha_18S_raw_hab <- as.data.frame(sapply(3:6, function(x) kruskal.test(alpha_div_18S_raw[,x],
                                                                       alpha_div_18S_raw[,1])))

# Rename columns with the proper variable names
kw_alpha_18S_raw_hab <- kw_alpha_18S_raw_hab %>%
  dplyr::rename(Observed = V1,
                Shannon = V2,
                Simpson = V3,
                Evenness = V4)

kw_alpha_18S_raw_hab <- t(kw_alpha_18S_raw_hab) # transpose
kw_alpha_18S_raw_hab <- as_tibble(kw_alpha_18S_raw_hab, rownames = "Metric") # adding rownames as a column
kw_alpha_18S_raw_hab <- kw_alpha_18S_raw_hab[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_18S_raw_hab) # checking object class




# KW analysis alpha_16S_nem on all metrics at once
# Remember, numerical variables are from columns 3-5. The test is by sitexhab, column 7
alpha_div_18S_raw_with_new_factor <- alpha_div_18S_raw

alpha_div_18S_raw_with_new_factor$Site_Habitat <- interaction(alpha_div_18S_raw_with_new_factor$Site, alpha_div_18S_raw_with_new_factor$Habitat)

kw_alpha_18S_raw_hab_site <- as.data.frame(sapply(3:6, function(x) kruskal.test(alpha_div_18S_raw_with_new_factor[,x],
                                                                                alpha_div_18S_raw_with_new_factor[,7])))

# Rename columns with the proper variable names
kw_alpha_18S_raw_hab_site <- kw_alpha_18S_raw_hab_site %>%
  dplyr::rename(Observed = V1,
                Shannon = V2,
                Simpson = V3,
                Evenness = V4)

kw_alpha_18S_raw_hab_site <- t(kw_alpha_18S_raw_hab_site) # transpose
kw_alpha_18S_raw_hab_site <- as_tibble(kw_alpha_18S_raw_hab_site, rownames = "Metric") # adding rownames as a column
kw_alpha_18S_raw_hab_site <- kw_alpha_18S_raw_hab_site[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_18S_raw_hab_site) # checking object class

# Save resulting table with fwrite to avoid any issues with characters
data.table::fwrite(kw_alpha_18S_raw, "kw_alpha_18S_raw.csv")
data.table::fwrite(kw_alpha_18S_raw_hab, "kw_alpha_18S_raw_hab.csv")
data.table::fwrite(kw_alpha_18S_raw_hab_site, "kw_alpha_18S_raw_hab_site.csv")


# Plot alpha diversity measures
# Change names

alpha_div_18S_raw$Site = factor(alpha_div_18S_raw$Site, levels=c("Campbell Cove","Westside Park","Mason's Marina"))

ad_18s_raw_site <- alpha_div_18S_raw %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Simpson", "Evenness")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Simpson", "Evenness"))) %>%
  ggplot(aes(x = Site, y = value)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  stat_compare_means() +
  stat_compare_means(comparisons = site_comparisons, p.adjust.methods = "BH", aes(label = ..p.signif..), size = 4, hide.ns = FALSE) +
  geom_jitter(aes(color = Site), height = 0, width = .2) +
  facet_nested(metric ~ Habitat, scales = "free") +
  scale_color_manual(values = alpha_color_site) +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_blank()) + # adjusts the title of x axis +
  scale_y_continuous(expand = c(0.1, 0.1)) +
  scale_x_discrete(
    labels = c("CC", "WP", "MM"),
    expand = c(0.2, 0.2),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) +  # adjusts the title of the legend
  ylab("") + # add the title on y axis
  xlab("Site") +  # add the title on x axis
  theme(legend.position="none")  # add the title on x axis

ad_18s_raw_site


ad_18s_raw_habitat <- alpha_div_18S_raw %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Simpson", "Evenness")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Simpson", "Evenness"))) %>%
  ggplot(aes(x = Habitat, y = value)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  stat_compare_means() +
  stat_compare_means(comparisons = habitat_comparisons, p.adjust.methods = "BH", aes(label =..p.signif..), size = 4, hide.ns = FALSE) +
  geom_jitter(aes(color = Habitat), height = 0, width = .2) +
  facet_nested(metric ~ Site, scales = "free") +
  scale_color_manual(values = alpha_color_habitat) +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_blank()) + # adjusts the title of x axis +
  scale_y_continuous(expand = c(0.1, 0.1)) +
  scale_x_discrete(
    labels = c("BS", "SG"),
    expand = c(0.2, 0.2),
    drop = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(strip.text = element_text(face = "bold", size =12), legend.title.align = 0.5) +  # adjusts the title of the legend
  ylab("") + # add the title on y axis
  xlab("Site") +  # add the title on x axis
  theme(legend.position="none")

ad_18s_raw_habitat
