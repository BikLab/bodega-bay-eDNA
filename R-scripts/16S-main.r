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
taxonomy_16s <- qiime2R::read_qza("16s/Updated QIIME2 Outputs/ca-16s-rep-sequences-taxonomy.qza")

table_16s <- qiime2R::read_qza("16s/Updated QIIME2 Outputs/16s-dada2-table-ca.qza")

metadata_16s <- read.csv("16s/Raw Data/16s_metadata.csv",
                         sep = ",", header = T, row = 1, quote = "")

tree_16s <- qiime2R::read_qza("16s/Updated QIIME2 Outputs/ca-rooted-16s-tree.qza")

# Extracted taxonomy data from taxonomy file
taxonomy_info_16s <- taxonomy_16s$data

# Parsed the taxonomy into 6 different taxonomic levels

parse_taxonomy_16s <- taxonomy_info_16s %>% separate(Taxon, sep =";", into = c("V1","V2","V3","V4","V5","V6", "V7"))

# Made Feature IDs the row names and removed the Feature.ID column 
rownames(parse_taxonomy_16s) <- parse_taxonomy_16s$Feature.ID 

parse_taxonomy_16s$Feature.ID <- NULL 

parse_taxonomy_16s$Consensus <- NULL 

# Removed any digits or characters from the beginning of taxon names
parse_taxonomy_16s[] <- lapply(parse_taxonomy_16s, function(x) gsub("D_\\d+__", "", x))

parse_taxonomy_16s[parse_taxonomy_16s==""] <- NA # Empty cells were classified as NA

parse_taxonomy_16s <- parse_taxonomy_16s %>% replace(is.na(.), "Unassigned") # Replaced all NAs with "Unassigned"

write.csv(parse_taxonomy_16s,"16s/Raw Data/taxonomy_ca.csv") # Exported parsed taxonomy as csv

# Transformed parsed taxonomy, table, and rooted tree into phyloseq objects
TAX_16s <- phyloseq::tax_table(as.matrix(parse_taxonomy_16s))

OTUMAT_16s <- table_16s$data

OTU_16s <- otu_table(OTUMAT_16s, taxa_are_rows = TRUE, replace_na(0))

otu_16s_export <- as.data.frame(OTU_16s)

write.csv(otu_16s_export,"16s_asv_table.csv") # Exported unfiltered (before decontam) ASV table as csv

TREE_16s <- tree_16s$data

# Creating phyloseq object
phylo_16s <- phyloseq(OTU_16s, sample_data(metadata_16s), TAX_16s, TREE_16s) # 9226 taxa and 100 samples

# Remove samples with low read counts and possible errors
to_remove <- c("16S-bodega-bay-fecal-experiment_CC.SG.3.1") # Possible error

to_remove2 <- c("16S-bodega-bay-fecal-experiment_MM.SG.1.3") # Low Read count

phylo_16s <- prune_samples(!(sample_names(phylo_16s) %in% to_remove), phylo_16s)

phylo_16s <- prune_samples(!(sample_names(phylo_16s) %in% to_remove2), phylo_16s)

# Remove samples with no reads
phylo_16s <- phylo_16s %>% 
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .)

# Subset to exclude Mitochondria, Chloroplasts, and Unassigned taxa at the domain level
phylo_16s <- subset_taxa(phylo_16s, V4 != "Chloroplast")

phylo_16s <- subset_taxa(phylo_16s, V5 != "Mitochondria")

phylo_16s <- subset_taxa(phylo_16s, V1 != "Unassigned")  # 9088 taxa and 97 samples

# Plot simple NMDS to view where samples, blanks, and controls fall
# Create new phyloseq object to provide transformation
phylo_16s_all_samples <- phylo_16s
phylo_16s_all_samples <- microbiome::transform(phylo_16s_all_samples, "compositional")

simple_ord_16s <- ordinate(phylo_16s_all_samples, "NMDS", "bray")
plot_ordination(phylo_16s_all_samples, simple_ord_16s, type="samples", color="Site", shape="Site")



###################### Using Decontam for filering out contaminants ######################
sample_data(phylo_16s)$is.neg <- sample_data(phylo_16s)$Habitat == "Control" | sample_data(phylo_16s)$Habitat == "Blank" # create a sample-variable for contaminants

phylo_contaminants <- isContaminant(phylo_16s, method = "prevalence", neg="is.neg", threshold=0.5, detailed = TRUE, normalize = TRUE) # detect contaminants based on control samples and their ASV prevalance

table(phylo_contaminants$contaminant) # check number of ASVs that are contaminents (75)

# Make phyloseq object of presence-absence in negative controls and true samples
phylo_contaminants.pa <- transform_sample_counts(phylo_16s, function(abund) 1 * (abund > 0)) # convert phyloseq table to presence-absence

ps.pa.neg <- prune_samples(sample_data(phylo_contaminants.pa)$Habitat == "Control" | sample_data(phylo_contaminants.pa)$Habitat == "Blank", phylo_contaminants.pa) # identify controls

ps.pa.pos <- prune_samples(sample_data(phylo_contaminants.pa)$Habitat != "Control" | sample_data(phylo_contaminants.pa)$Habitat != "Blank", phylo_contaminants.pa) # identify samples

df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg), contaminant=phylo_contaminants$contaminant) # convert into a dataframe

# Make phyloseq object of presence-absence in negative controls and true samples
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

phylo_16s <- prune_taxa(!phylo_contaminants$contaminant, phylo_16s) # remove ASVs identified as decontaminants from the dataset (9013 taxa and 97 samples)

clean_otu <- as.data.frame(otu_table(phylo_16s)) # Save ASV table as dataframe 

write.csv(clean_otu,"16s_clean_otu.csv") # Export ASV table as csv 



###################### Subsetting original phyloseq object for specific variables and normalizing ######################
# Subset any mocks, blanks, or controls that were left in the original phyloseq and normalize
phylo_16s <- phylo_16s %>% subset_samples(Habitat %in% c("Mason's Marina", "Campbell Cove", "West Park"))

phylo_normalized_16s <- phylo_normalized_16s <- transform_sample_counts(phylo_16s, function(x) x/sum(x)) 

# Subset phyloseq by habitat and normalize
phylo_16s <- phylo_16s %>% subset_samples(Site %in% c("Campbell Cove", "Westside Park", "Mason's Marina"))

phylo_normalized_16s <- microbiome::transform(phylo_16s, "compositional") 

# Subset phyloseq by site and normalize
CC_phylo_16s <- phylo_16s %>% subset_samples(Site %in% ("Campbell Cove"))
CC_phylo_16s_normalized <- microbiome::transform(CC_phylo_16s, "compositional") 

WP_phylo_16s <- phylo_16s %>% subset_samples(Site %in% ("Westside Park"))
WP_phylo_16s_normalized <- microbiome::transform(WP_phylo_16s, "compositional")

MM_phylo_16s <- phylo_16s %>% subset_samples(Site %in% ("Mason's Marina"))
MM_phylo_16s_normalized <- microbiome::transform(MM_phylo_16s, "compositional")



###################### Plotting for Relative Abundance ######################
# Functions to collapse a certain number of taxa into category others 

# Top 10 most abundant taxa + others
merge_top10 <- function(ps_object, top=9){
  transformed <- transform_sample_counts(ps_object, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:7] <- "Others"} # 1:7 if there are species level
  }
  return(merged)
}

# Top 20 most abundant taxa + others
merge_top20 <- function(ps_object, top=19){
  transformed <- transform_sample_counts(ps_object, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:7] <- "Others"} # 1:7 if there are species level
  }
  return(merged)
}

# Generate colors for taxa (top 10 & 20)
colors_top10 <- c("#9a6324", "#46f0f0", "#aaffc3", "#33A02C", "#f032e6", "#bcf60c", "#f58231", "#6A3D9A", "#fffac8", "gray")

colors_top20 <- c("#dd3497", "#ae017e","#7a0177","#3690c0", "#74a9cf", "#000075", "#a6bddb", "#d0d1e6", "#014636", "#016c59", 
"#02818a", "#41b6c4", "#7fcdbb","#c7e9b4","#e0f3db", "#ccece6",  "#f768a1", "#fa9fb5", "#fcc5c0", "gray")

# Generating phyloseq objects for top 10 taxa
temp_phylo_16s <- phylo_16s # Use different variable name to avoid overwriting original phyloseq object

# Use tax_fix function to assign taxa names to lower ranks that are unknown
fix_phylo_16s <- tax_fix(temp_phylo_16s, unknowns = c("Unassigned", "uncultured", "Unkown", "uncultured bacterium", "Unknown Family")) 

# Agglomerate taxa down to V4 rank (Order level) and V5 rank (Family level)
glom_16s <- tax_glom(fix_phylo_16s, taxrank = "V4")

glom_16s_fam <- tax_glom(fix_phylo_16s, taxrank = "V5")

# Top 10 at order level
# Run function on phyloseq object and make into data frame
phy_16_top10_V4 <- merge_top10(glom_16s, top=9)

phy_16_top10_V4_df <- psmelt(phy_16_top10_V4)

# Add common factors to use for plotting
phy_16_top10_agr = aggregate(Abundance~Sample+Site+Habitat+V4, data=phy_16_top10_V4_df, FUN=mean) 

# Reorder taxa to put "Others" as the final taxa of the Order list
unique(phy_16_top10_agr$V4)
phy_16_top10_agr$V4 <- factor(phy_16_top10_agr$V4,
                              levels = c("Campylobacterales", "Cellvibrionales","Desulfobacterales","Desulfobulbales","Flavobacteriales", 
                                         "Gammaproteobacteria V3", "Pirellulales","Rhizobiales", "Rhodobacterales", "Others"))

# Reorder Site levels
phy_16_top10_agr$Site = factor(phy_16_top10_agr$Site, levels=c("Campbell Cove", "Westside Park", "Mason's Marina"))

# Plot by site and habitat - Order level top 10
taxonomy_bar_16s_top10_order <- ggplot(phy_16_top10_agr, aes(x = Sample, y = Abundance, fill = V4)) +
  facet_nested(. ~ Site+Habitat, scales = "free") + # facet by Site and Habitat
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 1, paste(round(Abundance*100, digits = 0), "%"), "")), size = 2, position = position_stack(vjust = 0.5)) + # adds text labels to the bars
  scale_fill_manual(values = colors_top10, name = "Order") + # Uses custom colors for the taxa and give title to the legend
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, vjust = 0.5, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) + # Format facet grid title 
  scale_x_discrete(label = function(x) stringr::str_replace(x, "16S-bodega-bay-fecal-experiment_","") # remove the prefix from the sample names
  ) +

taxonomy_bar_16s_top10_order



# Top 10 at family level
# Run function on phyloseq object and make into data frame
phy_16_top10_V5 <- merge_top10(glom_16s_fam, top=9)

phy_16_top10_V5_df <- psmelt(phy_16_top10_V5)

# Add common factors to use for plotting
phy_16_top10_agr_fam = aggregate(Abundance~Sample+Site+Habitat+V5, data=phy_16_top10_V5_df, FUN=mean) 

# Reorder taxa to put "Others" as the final taxa of the Family list
unique(phy_16_top10_agr_fam$V5)
phy_16_top10_agr_fam$V5 <- factor(phy_16_top10_agr_fam$V5,
                                  levels = c("Anaerolineaceae", "Desulfocapsaceae","Desulfosarcinaceae","Flavobacteriaceae","Gammaproteobacteria V3",
                                             "Halieaceae", "Pirellulaceae","Rhodobacteraceae","Sulfurovaceae","Others"))

# Reorder Site levels
phy_16_top10_agr_fam$Site = factor(phy_16_top10_agr_fam$Site, levels=c("Campbell Cove","Westside Park","Mason's Marina"))

# Plot by site and habitat - Family level top 10
taxonomy_bar_16s_top10_family <- ggplot(phy_16_top10_agr_fam, aes(x = Sample, y = Abundance, fill = V5)) +
  facet_nested(. ~ Site+Habitat, scales = "free") + # facet by Site and Habitat
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 1, paste(round(Abundance*100, digits = 0), "%"), "")), size = 2, position = position_stack(vjust = 0.5)) + # adds text labels to the bars
  scale_fill_manual(values = colors_top10, name = "Family") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, vjust = 0.5, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) + # Format facet grid title 
  scale_x_discrete(label = function(x) stringr::str_replace(x, "16S-bodega-bay-fecal-experiment_","")) # remove the prefix from the sample names

taxonomy_bar_16s_top10_family



# Top 20 at order level
# Run function on phyloseq object and make into data frame
phy_16_top20_V4 <- merge_top20(glom_16s, top=19)

phy_16_top20_V4_df <- psmelt(phy_16_top20_V4)

# Add common factors to use for plotting
phy_16_top20_agr = aggregate(Abundance~Sample+Site+Habitat+V4, data=phy_16_top20_V4_df, FUN=mean) 

# Put "Others" to the final of the Order list - top 20
unique(phy_16_top20_agr$V4)
phy_16_top20_agr$V4 <- factor(phy_16_top20_agr$V4,
                                   levels = c("Actinomarinales", "Anaerolineales","B2M28","Bacteroidales","Campylobacterales", "Cellvibrionales",
                                              "Chromatiales","Desulfobacterales","Desulfobulbales", "Ectothiorhodospirales", "Flavobacteriales", 
                                              "Gammaproteobacteria Incertae Sedis", "Gammaproteobacteria V3", "Pirellulales","Rhizobiales", 
                                              "Rhodobacterales","Steroidobacterales","Thiotrichales", "Verrucomicrobiales", "Others"))

# Reorder Site levels
phy_16_top20_agr$Site = factor(phy_16_top20_agr$Site, levels=c("Campbell Cove","Westside Park","Mason's Marina"))

# Plot by site - Order level top 10
taxonomy_bar_16s_order <- ggplot(phy_16_top20_agr, aes(x = Sample, y = Abundance, fill = V4)) +
  facet_nested(. ~ Site+Habitat, scales = "free") + # facet by Site and Habitat
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 2, position = position_stack(vjust = 0.5)) +  # adds text labels to the bars 
  scale_fill_manual(values = colors_top20_3, name = "Order") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, vjust = 0.5, face = "bold")) +# adjusts text of y axis
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
  scale_x_discrete(label = function(x) stringr::str_replace(x, "16S-bodega-bay-fecal-experiment_","")) #  remove the prefix from the sample names

taxonomy_bar_16s_order



# Top 20 at family level
# Run function on phyloseq object and make into data frame
phy_16_top20_V5 <- merge_top20(glom_16s_fam, top=19)

phy_16_top20_V5_df <- psmelt(phy_16_top20_V5)

# Add common factors to use for plotting
phy_16_top20_agr_fam = aggregate(Abundance~Sample+Site+Habitat+V5, data=phy_16_top20_V5_df, FUN=mean) 
unique(phy_16_top20_agr_fam$V5)

# Put "Others" to the final of the Family list - top 20
phy_16_top20_agr_fam$V5 <- factor(phy_16_top20_agr_fam$V5,
                              levels = c("Actinomarinales V4", "Anaerolineaceae","B2M28 V4","Bacteroidetes BD2-2","Desulfobulbaceae", "Desulfocapsaceae",
                                         "Desulfosarcinaceae","Flavobacteriaceae","Gammaproteobacteria Incertae Sedis V4", "Gammaproteobacteria V3", "Halieaceae", 
                                         "Methyloligellaceae", "Milano-WF1B-44 V4", "Pirellulaceae","Rhodobacteraceae", 
                                         "Sedimenticolaceae","Sulfurovaceae","Thiotrichaceae", "Woeseiaceae", "Others"))

# Reorder Site levels
phy_16_top20_agr_fam$Site = factor(phy_16_top20_agr_fam$Site, levels=c("Campbell Cove","Westside Park","Mason's Marina"))

# Plot by site - Family level top 20
taxonomy_bar_16s_family <- ggplot(phy_16_top20_agr_fam, aes(x = Sample, y = Abundance, fill = V5)) +
  facet_nested(. ~ Site+Habitat, scales = "free") + # facet by Site and Habitat
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 2, position = position_stack(vjust = 0.5)) + # adds text labels to the bars
  scale_fill_manual(values = colors_top20_3, name = "Family") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, vjust = 0.5, face = "bold")) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_text(face = "bold", size = 12)) + # adjusts the title of x axis
  scale_y_continuous(labels=scales::percent, expand = c(0.0, 0.0)) + # plot as % and removes the internal margins
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white")) + # removes the gridlines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # adjusts the title of the legend
  ylab("Relative Abundance") + # add the title on y axis
  xlab("Sample") + # add the title on x axis
  theme(strip.background =element_rect(
    color = "black",
    fill = "white",
    linewidth = 1,
    linetype = "solid"),
    strip.text = element_text(
      size = 12, color = "black", face = "bold")) + # Format facet grid title 
  scale_x_discrete(label = function(x) stringr::str_replace(x, "16S-bodega-bay-fecal-experiment_","")) # remove the prefix from the sample names

taxonomy_bar_16s_family



###################### Beta Diversity Ordination ######################
# Create ordination object using PCoA and Bray-Curtis dissimilarity
phylo_16s_ord <- ordinate(phylo_normalized_16s, "PCoA", "bray") 

# Plot the ordination object with samples colored by Site and shaped by Habitat
phylo_16s_pcoa <- plot_ordination(phylo_normalized_16s, phylo_16s_ord, type="samples", color="Site", shape="Habitat", title = "16S rRNA - Microbes") +
  geom_point(size = 3) +
  scale_color_manual(values = c("#C8AB83","#2E5266","#E43F6F"),breaks = c("Campbell Cove", "Westside Park", "Mason's Marina")) + #  custom colors for Site
  scale_shape_manual(values = c(1,19), name = "Habitat", breaks = c("Sea Grass", "Bare Sediment")) + # custom shapes for Habitat
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 14, color = "black", face = "bold")) + # adjusts text of y axis
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,size = 12, color = "black")) + # adjusts text of y axis
  theme(legend.title = element_text(face = "bold")) + # adjusts the title of the legend
  theme(axis.line = element_line(color = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # removes the gridlines
  theme(plot.title = element_text(color = "black", face = "bold")) + # adjusts the title of the plot
  theme(legend.position = "none") # removes the legend for when needing to combine with other plots

phylo_16s_pcoa$layers <- phylo_16s_pcoa$layers[-1] # Removes the points layer to avoid duplication

phylo_16s_pcoa

# PERMANOVA for 16S Samples
# Calculate bray curtis distance matrix
bray_16s <- phyloseq::distance(phylo_normalized_16s, method = "bray")

# Make a data frame from the sample_data
sampledf_16s <- data.frame(sample_data(phylo_16s))

# Adonis test
adonis_results_16s_site <- adonis2(bray_16s ~ Site, data = sampledf_16s)
adonis_results_16s_site # Significance found at Pr(>F) = .001 (***)

# Homogeneity of dispersion test
beta_16s_results_site <- betadisper(bray_16s, sampledf_16s$Site)
permutest(beta_16s_results_site) # Significance found at Pr(>F) = .001 (***)

# Pairwise Wilcoxon test for 16S samples by Site
pairwise_16s_site <- pairwise.wilcox.test(bray_16s, sampledf_16s$Site, 
                                              p.adjust.method = "bonferroni")


                                              
###################### Alpha Diversity ######################
# Calculate alpha-diversity measures (For plot purposes only!)
# This can be done using the different phyloseq alpha diversity measures
# You will get a Warning message for each index since there is no singletons on the dataset

# Sample 16S-bodega-bay-fecal-experiment_CC.SG.3.2 was removed as an outlier due to high read count
to_remove3 <- c("16S-bodega-bay-fecal-experiment_CC.SG.3.2")

phylo_16s_ad <- prune_samples(!(sample_names(phylo_16s) %in% to_remove3), phylo_16s)

# Calculate alpha diversity measures for 16S data
alpha_div_16S <- data.frame(
  "Observed" = phyloseq::estimate_richness(phylo_16s_ad, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(phylo_16s_ad, measures = "Shannon"),
  "InvSimpson" = phyloseq::estimate_richness(phylo_16s_ad, measures = "InvSimpson"),
  "Site" = phyloseq::sample_data(phylo_16s_ad)$Site,
  "Habitat" = phyloseq::sample_data(phylo_16s_ad)$Habitat)

# Add Evenness measure to the alpha_div_16S dataframe
alpha_div_16S$Evenness <- alpha_div_16S$Shannon/log(alpha_div_16S$Observed)

head(alpha_div_16S)

# Rename variable InvSimpson to Simpson
# The function rename & %>% works on dplyr. make sure it is loaded.
alpha_div_16S <- alpha_div_16S %>%
  dplyr::rename(Simpson = InvSimpson)
head(alpha_div_16S)

# Reorder dataframe, first categorical then numerical variables  
alpha_div_16S <- alpha_div_16S[, c( 5, 4, 1, 2, 3, 6)]
head(alpha_div_16S)

#Summarize alpha diversity measures by site
summary_alpha_16S_site <- alpha_div_16S %>%
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

write_csv(summary_alpha_16S_site, "16s/Tables:Results/Alpha Diversity/alpha_div_results_16s_site.csv") # save summary table

#Summarize alpha diversity measures by habitat
summary_alpha_16S_habitat <- alpha_div_16S %>%
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

write_csv(summary_alpha_16S_habitat, "16s/Tables:Results/Alpha Diversity/alpha_div_results_16s_habitat.csv") #  save summary table

#Summarize alpha diversity measures by habitat within site
summary_alpha_16S_site_habitat <- alpha_div_16S %>%
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

write_csv(summary_alpha_16S_site_habitat, "16s/Tables:Results/Alpha Diversity/summary_alpha_16S_site_habitat.csv") # save summary table


# KW analysis alpha_div_16S on all metrics at once
# Remember, numerical variables are from columns 3-6. The test is by site, column 2
kw_alpha_16S_site <- as.data.frame(sapply(3:6, function(x) kruskal.test(alpha_div_16S[,x], alpha_div_16S[,2])))

# Rename columns with the proper variable names
kw_alpha_16S_site <- kw_alpha_16S_site %>%
  dplyr::rename(Observed = V1,
         Shannon = V2,
         Simpson = V3,
         Evenness = V4)

kw_alpha_16S_site <- t(kw_alpha_16S_site) # transpose
kw_alpha_16S_site <- as_tibble(kw_alpha_16S_site, rownames = "Metric") # adding rownames as a column
kw_alpha_16S_site <- kw_alpha_16S_site[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_16S_site) # checking object class


# Save resulting table with fwrite to avoid any issues with characters
data.table::fwrite(kw_alpha_16S_site, "16s/Tables:Results/Alpha Diversity/kw_ad_16s_site.csv")

# KW analysis alpha_div_16S on all metrics at once
# Remember, numerical variables are from columns 3-5. The test is by habitat, column 1
kw_alpha_16S_habitat <- as.data.frame(sapply(3:6, function(x) kruskal.test(alpha_div_16S[,x], alpha_div_16S[,1])))

# Rename columns with the proper variable names
kw_alpha_16S_habitat <- kw_alpha_16S_habitat %>%
  dplyr::rename(Observed = V1,
         Shannon = V2,
         Simpson = V3,
         Evenness = V4)

kw_alpha_16S_habitat <- t(kw_alpha_16S_habitat) # transpose
kw_alpha_16S_habitat <- as_tibble(kw_alpha_16S_habitat, rownames = "Metric") # adding rownames as a column
kw_alpha_16S_habitat <- kw_alpha_16S_habitat[, -(5:6)] # removing columns 5 and 6
class(kw_alpha_16S_habitat) # checking object class


# Save resulting table with fwrite to avoid any issues with characters
data.table::fwrite(kw_alpha_16S_habitat, "16s/Tables:Results/Alpha Diversity/kw_ad_16s_habitat.csv")

# Plot alpha diversity measures
alpha_div_16S$Site = factor(alpha_div_16S$Site, levels=c("Campbell Cove","Westside Park","Mason's Marina"))

# Create a list of comparisons for the sites
# This will be used for the statistical comparisons in the plots
site_comparisons <- list( c("Campbell Cove", "Westside Park"), c("Westside Park", "Mason's Marina"), c("Campbell Cove", "Mason's Marina"))

# Create custom colors for the sites
alpha_color_site <- c("#C8AB83","#2E5266","#E43F6F")

ad_16s_site <- alpha_div_16S %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Simpson", "Evenness")) %>% # gather the metrics into long format
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Simpson", "Evenness"))) %>% # reorder the metrics
  ggplot(aes(x = Site, y = value)) +  
  geom_boxplot(outlier.color = NA, width = 0.5) + # creates the boxplot
  stat_compare_means(comparisons = site_comparisons, p.adjust.methods = "BH", aes(label = ..p.signif..), size = 4, hide.ns = FALSE) + # adds the statistical comparisons
  geom_jitter(aes(color = Site), height = 0, width = .2) + # adds the points to the boxplot
  facet_nested(metric ~ Habitat, scales = "free") + # facets by metric and habitat
  scale_color_manual(values = alpha_color_site) + # uses custom colors for the sites
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
  theme(legend.position="none")  # add/removes the title on x axis

ad_16s_site



# Create a list of comparisons for the habitats
alpha_color_habitat <- c("saddlebrown","#00A572")

# Create custom colors for the habitats
habitat_comparisons <- list(c("Bare Sediment", "Sea Grass"))

ad_16s_habitat <- alpha_div_16S %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Simpson", "Evenness")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Simpson", "Evenness"))) %>%
  ggplot(aes(x = Habitat, y = value)) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  stat_compare_means(method = "wilcox.test", comparisons = habitat_comparisons, p.adjust.methods = "BH", aes(label =..p.signif..), size = 4, hide.ns = FALSE) +
  geom_jitter(aes(color = Habitat), height = 0, width = .2) +
  facet_nested(~ metric, scales = "free") +
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

ad_16s_habitat