####################################### Using Phyloseq #########################################
# Load table and make adjustments 
ec_number <- read.table("Picrust/Picrust After/pred_metagenome_unstrat_descrip.tsv", header=T, 
                        sep="\t", stringsAsFactors=F, quote = "", check.names=F, comment.char="")

ec_number <- ec_number[, -(2)]

rownames(ec_number) <- ec_number$`function`

ec_number$`function` <- NULL 

# Round all numeric columns to 2 decimal places
ec_number_rounded <- ec_number %>% 
  mutate(across(where(is.numeric), ~ round(., digits = 0)))

# Create table with gene function ids and description only
ec_number_descrip <- read.table("Picrust/Picrust After/pred_metagenome_unstrat_descrip.tsv", header=T, 
                                sep="\t", stringsAsFactors=F, quote = "", check.names=F, comment.char="")

ec_number_descrip <- ec_number_descrip[, -(3:102)]

ec_number_descrip <- ec_number_descrip %>%
  rename(EC_number = `function`)

# Load metadata files for each experiment
metadata <- read.csv("16s/Raw Data/16s_metadata.csv", sep = ",", header = T, row = 1, quote = "")

# Creating phyloseq object
metadata_ec <- ec_number_descrip
metadata_sample <- metadata
ec_predicted <- ec_number_rounded


rownames(metadata_ec) <- metadata_ec$EC_number
metadata_ec$EC_number <- NULL 


ec_phy <- otu_table(ec_predicted, taxa_are_rows = TRUE) # notes taxa (AVSs) are as rows
tax_phy <- tax_table(as.matrix(metadata_ec))
samples <- sample_data(metadata_sample)

write.csv(metadata_sample, "16S-EC-metadata.csv")
write.csv(ec_phy, "16S-EC-table.csv")
write.csv(metadata_ec, "16S-EC-table-description.csv")


picrust_phyloseq <- phyloseq(ec_phy, tax_phy, samples) # 3046 taxa and 100 samples

####################################### Making edits to phyloseq object #######################################

combined_ec_asv_table <- read.table("Picrust/Picrust After/combined_EC_predicted.tsv", header=T, 
                                    sep="\t", stringsAsFactors=F, quote = "", check.names=F, comment.char="")

names(combined_ec_asv_table)[1] <- "OTU"

rownames(combined_ec_asv_table) <- combined_ec_asv_table$OTU
combined_ec_asv_table <- combined_ec_asv_table[, -1]

combined_ec_asv_matrix <- t(combined_ec_asv_table)

combined_ec_asv_df <- as.data.frame(combined_ec_asv_matrix)

# Attempt at new phylo object with intermediate phylo
phylo_16s_fixed <- tax_fix(phylo_16s, unknowns = c("Unassigned", "uncultured", "Unkown", "uncultured bacterium", 
                                                   "Unknown Family", "uncultured organism", "uncultured gamma proteobacterium", 
                                                   "uncultured delta proteobacterium", "unknown archaeon", "uncultured archaeon",
                                                   "uncultured crenarchaeote", "uncultured euryarchaeote"))

intermediate_asv_table <- otu_table(combined_ec_asv_table, taxa_are_rows = TRUE) # notes taxa (AVSs) are as rows
intermediate_tax <- tax_table(phylo_16s_fixed)
intermediate_metadata <- sample_data(metadata_ec)

intermediate_phyloseq <- phyloseq(intermediate_asv_table, intermediate_tax, intermediate_metadata)

# create individual extractions
otu_mat_intermediate <- as(otu_table(intermediate_phyloseq), "matrix")
tax_mat_intermediate <- as(tax_table(intermediate_phyloseq), "matrix")
# make edit to column names
write.csv(metadata_ec, "intermediate_metadata.csv")
ec_meta <- read.csv("intermediate_metadata.csv", row.names = 1)



##########################################Find EC numbers with a specific keyword #######################################
matching_ecs <- rownames(metadata_ec)[grepl("sulf", metadata_ec$description, ignore.case = TRUE)]

# Initialize a list to hold results for each EC
result_list <- lapply(matching_ecs, function(ec) {
  # Find which OTUs are associated with this EC
  otus_for_ec <- rownames(otu_mat_intermediate)[otu_mat_intermediate[, ec] > 0]
  # Get the taxonomy assignments for these OTUs
  tax_for_ec <- tax_mat_intermediate[otus_for_ec, , drop = FALSE]
  # Return as a list with OTU IDs and taxonomy table
  list(otus = otus_for_ec, taxonomy = tax_for_ec)
})
names(result_list) <- matching_ecs


summarize_genus_row <- function(tax_tab) {
  # Get the Genus column; adjust "Genus" if your column is named differently
  as.character(tax_tab[, "V6"])
}


summary_table <- do.call(rbind, lapply(names(result_list), function(ec) {
  row <- result_list[[ec]]
  if (NROW(row$taxonomy) == 0) {
    genus_str <- ""
  } else {
    genus_str <- paste(unique(na.omit(summarize_genus_row(row$taxonomy))), collapse = ", ")
  }
  data.frame(
    EC_Number = ec,
    Description = ec_meta[ec, "description"],
    Contributing_Genera = genus_str,
    stringsAsFactors = FALSE
  )
}))

write.csv(summary_table, "summary_table_with_bacteria_genus.csv")

# Expand all genera into a single vector
all_genera <- unlist(strsplit(summary_table$Contributing_Genera, ",\\s*"))

# Count how many times each genus appears
genus_counts <- as.data.frame(table(all_genera), stringsAsFactors = FALSE)

# Sort in descending order of frequency
top_genera <- genus_counts %>% arrange(desc(Freq))

write.csv(top_genera, "top_genera_for_sulfur_ec.csv")

# Get the top 20 most frequent genera
top_20_genera <- head(top_genera, 20)

top_30_genera <- head(top_genera, 30)

top_40_genera <- head(top_genera, 40)

top_50_genera <- head(top_genera, 50)

top_100_genera <- head(top_genera, 100)


picrust_phyloseq <- subset_samples(picrust_phyloseq, 
                                   Site %in% c("Campbell Cove", "Westside Park", "Mason's Marina")) # 3046 taxa and 84 samples

picrust_phyloseq <- prune_samples(!(sample_names(picrust_phyloseq) %in% to_remove), picrust_phyloseq)

picrust_phyloseq <- prune_samples(!(sample_names(picrust_phyloseq) %in% to_remove2), picrust_phyloseq)

# Remove samples with low read counts and possible errors
picrust_phyloseq_prune <- picrust_phyloseq %>% 
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .) # 2974 taxa and 82 samples

picrust_phyloseq_relab <- microbiome::transform(picrust_phyloseq_prune, "compositional")

picrust_df <- psmelt(picrust_phyloseq_prune)

######################################## Ordinations ################################
picrust_phyloseq_pcoa <- ordinate(picrust_phyloseq_relab, "PCoA", "bray") 
picrust_ord_plot <- phyloseq::plot_ordination(picrust_phyloseq_relab, picrust_phyloseq_pcoa, type="samples", color="Site", 
                                    shape="Habitat", title = "Picrust Predictive Functions") + 
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
  scale_y_reverse() +
  scale_x_reverse()

picrust_ord_plot$layers <- picrust_ord_plot$layers[-1]

picrust_ord_plot

# Stat test
# PERMANOVA for All 16s Samples
# Calculate bray curtis distance matrix
bray_ec <- phyloseq::distance(picrust_phyloseq_relab, method = "bray")

# make a data frame from the sample_data
sampledf_ec <- data.frame(sample_data(picrust_phyloseq_relab))

# Adonis test
adonis_results_ec <- adonis2(bray_ec ~ Site, data = sampledf_ec)
adonis_results_ec # Significance found at Pr(>F) = .001 (***)

# Homogeneity of dispersion test
beta_results_ec <- betadisper(bray_ec, sampledf_ec$Site)
  permutest(beta_results_ec) # Significance found at Pr(>F) = .001 (***)


picrust_phyloseq_cc_normalized <- microbiome::transform(picrust_phyloseq_cc, "compositional")

picrust_cc_pcoa <- ordinate(picrust_phyloseq_cc_normalized, "PCoA", "bray") 

picrust_cc_plot <- plot_ordination(picrust_phyloseq_cc_normalized, picrust_cc_pcoa, type="samples", color="Habitat", title = "16S rRNA Campbell Cove - Microbes") + 
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

picrust_cc_plot$layers <- picrust_cc_plot$layers[-1]

picrust_cc_plot




picrust_phyloseq_mm_normalized <- microbiome::transform(picrust_phyloseq_mm, "compositional")
picrust_mm_pcoa <- ordinate(picrust_phyloseq_mm_normalized, "PCoA", "bray") 
picrust_mm_plot <- plot_ordination(picrust_phyloseq_mm_normalized, picrust_mm_pcoa, type="samples", color="Habitat", title = "16S rRNA Mason's Marina - Microbes") + 
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

picrust_mm_plot$layers <- picrust_mm_plot$layers[-1]

picrust_mm_plot


picrust_phyloseq_wp_normalized <- microbiome::transform(picrust_phyloseq_wp, "compositional")
picrust_wp_pcoa <- ordinate(picrust_phyloseq_wp_normalized, "PCoA", "bray") 
picrust_wp_plot <- plot_ordination(picrust_phyloseq_wp_normalized, picrust_wp_pcoa, type="samples", color="Habitat", title = "16S rRNA Westside Park - Microbes") + 
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

picrust_wp_plot$layers <- picrust_wp_plot$layers[-1]

picrust_wp_plot


################################################# Bar Plot ############################
top20functions <- names(sort(taxa_sums(picrust_phyloseq_relab), TRUE)[1:20])

#subset phyloseq object to only selected taxa
top20functions_phy <- prune_taxa(top20functions, picrust_phyloseq_relab_prune) 

top20functions_phy_df <- psmelt(top20functions_phy)

# Reorder Site levels
top20functions_phy_df$Site = factor(top20functions_phy_df$Site, levels=c("Campbell Cove","Westside Park","Mason's Marina"))


picrust_bar_plot <- ggplot(top20functions_phy_df, aes(x = Sample, y = Abundance, fill = description)) +
  facet_nested(. ~ Site+Habitat, scales = "free",
               labeller = labeller(Site = site.labs, Habitat = habitat.labs)) +
  geom_bar(stat = "identity", width = 0.95) + # adds to 100%
  geom_text(aes(label = ifelse(round(Abundance*100) >= 5, paste(round(Abundance*100, digits = 0), "%"), "")), size = 2, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors_top20_2, name = "Description") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face = "bold")) + # adjusts text of y axis
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + # adjusts text of x axis
  theme(axis.title.y = element_text(face = "bold", size = 12)) +  # adjusts the title of y axis
  theme(axis.title.x = element_blank()) + # adjusts the title of x axis
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
      size = 10, color = "black", face = "bold")) + # Format facet grid title 
  scale_x_discrete(label = function(x) stringr::str_replace(x, "16S-bodega-bay-fecal-experiment_","")) 

picrust_bar_plot



####################################### Run ALEDx2 for sites ###################################
aldex2_ec <- ALDEx2::aldex(data.frame(phyloseq::otu_table(picrust_phyloseq)),
                            phyloseq::sample_data(picrust_phyloseq)$Site,
                            mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
# Create data frame and format
ec_otu_table <- data.frame(phyloseq::otu_table(picrust_phyloseq))
ec_otu_table <- rownames_to_column(ec_otu_table, var = "OTU")

write.csv(aldex2_ec, "Picrust/Aldex2/aldex2_ec.csv")

# Import results CSV and format
aldex2_ec_result <- read.csv("Picrust/Aldex2/aldex2_ec.csv")
colnames(aldex2_ec_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

# Create data frame for taxonomy 
aldex_taxa_info_ec <- data.frame(tax_table(picrust_phyloseq))
aldex_taxa_info_ec <- aldex_taxa_info_ec %>%
  rownames_to_column(var = "OTU")

sample_tab_ec <- data.frame(sample_data(picrust_phyloseq))

# Filter aldex2 results by sig kw.ep and join the taxanomic information
sig_aldex2_ec_result <- aldex2_ec_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)

sig_aldex2_ec_result <- left_join(sig_aldex2_ec_result, aldex_taxa_info_ec)

#ec_number_aldex_sig_top20 <- sig_aldex2_ec_result %>% top_n(-20, kw.ep)

#ec_number_aldex_sig_top20 <- left_join(ec_number_aldex_sig_top20, aldex_taxa_info_ec)

write.csv(sig_aldex2_ec_result, "Picrust/Picrust After/sig_aldex2_ec_result.csv")


# Made changes to taxonomy names if needed 
#sig_aldex2_ec_result <- read.csv("Picrust/Aldex2/sig_aldex2_ec_result.csv")

sig_aldex2_ec_result <- sig_aldex2_ec_result[grepl("sulf", sig_aldex2_ec_result$description, ignore.case = TRUE), ]

ec_number_aldex_sig_top20 <- sig_aldex2_ec_result %>% top_n(-25, kw.ep)

# sig_aldex2_ec_result <- sig_aldex2_ec_result[-4,]


# Create clr objects by using OTU ids, original otu table, and all significant OTUs
sig_aldex2_ec_result_count <- left_join(ec_number_aldex_sig_top20, ec_otu_table)
#sig_aldex2_ec_result_count <- sig_aldex2_ec_result_count[, -1]
write.csv(sig_aldex2_ec_result_count, "Picrust/Aldex2/sig_aldex2_ec_result_count_sulfur.csv")

# Dropped taxonomic rank and formatted columns 
clr_ec <- sig_aldex2_ec_result_count[, -(2:4)] 
rownames(clr_ec) <- clr_ec$OTU
clr_ec <- clr_ec[, -1]

# Adjusting zeros on the matrix and applying log transformation
shsk_ec_czm <- cmultRepl(t(clr_ec),  label=0, method="CZM")
shsk_ec_czm_tv <- t(apply(shsk_ec_czm, 1, function(x){log(x) - mean(log(x))}))
shsk_ec_czm <- (apply(clr_ec, 1, function(x){log(x+1) - mean(log(x+1))}))


# Combine the heatmap and the annotation
Z.Score.ec <- scale(t(shsk_ec_czm))

# Organize total abundance and taxa name
# Taxa name
Z.Score.ec_name <- as.data.frame(Z.Score.ec)
str(Z.Score.ec_name)
Z.Score.ec_name <- rownames_to_column(Z.Score.ec_name, var = "OTU")
Z.Score.ec_name <- left_join(Z.Score.ec_name, aldex_taxa_info_ec)
head(Z.Score.ec_name)

# Taxa abundance 
ec_otu_table_total <- as.data.frame(ec_otu_table)
ec_otu_table_total$Total <- rowSums(ec_otu_table[, -1])
head(ec_otu_table_total)
ec_otu_table_total <- ec_otu_table_total[, -(2:83)] 

Z.Score.ec_count_total <- left_join(Z.Score.ec_name, ec_otu_table_total)
head(Z.Score.ec_count_total)

# Defining color scheme for row annotations
abundance_col_fun_ec = colorRamp2(c(0, 200000, 500000, 1000000),
                                   c("#c7eae5",
                                     "#80cdc1",
                                     "#35978f",
                                     "#01665e"))

ha_right_ec = rowAnnotation(
  Abundance = Z.Score.ec_count_total$Total, border = FALSE, col = list(Abundance = abundance_col_fun_ec))
row_labels_ec = Z.Score.ec_count_total$description


ha_ec = HeatmapAnnotation(
  Habitat = as.vector(sample_tab_ec$Habitat),
  Site = as.vector(sample_tab_ec$Site),
  col = list(
    Habitat = c("Bare Sediment" = "saddlebrown", "Sea Grass" = "#00A572"),
    Site = c("Mason's Marina" = "#E43F6F", "Westside Park" = "#2E5266",
             "Campbell Cove" = "#C8AB83")
  ),
  annotation_legend_param = list(
    Habitat = list(
      title = "Habitat",
      at = c("Bare Sediment", "Sea Grass"),
      labels = c ("BS", "SG")
    ),
    Site = list(
      title = "Site",
      at = c("Campbell Cove", "Westside Park", "Mason's Marina"),
      labels = c ("CC", "WP", "MM")
    )
  ))

# Plot heatmap at the EC level
hm_ec <- Heatmap(Z.Score.ec, name = "Z-score, CLR", col = col_matrix,
                  column_title  = "Sulfur EC Numbers", 
                  column_title_gp = gpar(fontface = "bold", fontsize = 14),
                  column_split = as.vector(as.vector(sample_tab_ec$Site)),
                  column_order = order(as.numeric(gsub("column", "", colnames(Z.Score.ec)))),
                  border = TRUE,
                  top_annotation = ha_ec,
                  right_annotation = ha_right_ec,
                  row_title = "EC",
                  row_labels = row_labels_ec,
                  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                  row_names_gp = gpar(fontsize = 6),
                  column_names_gp = gpar(fontsize = 6),
                  #row_order = order(row_labels_ec),
                  rect_gp = gpar(col = "white", lwd = 1),
                  show_column_names = FALSE,
                  show_heatmap_legend = TRUE)


hm_ec


 ####################################### Run ALEDx2 for Campbell Cove only ###################################


picrust_phyloseq_cc <- subset_samples(picrust_phyloseq, Site == "Campbell Cove")

d1 <- psmelt(picrust_phyloseq_cc)

aldex2_ec_cc <- ALDEx2::aldex(data.frame(phyloseq::otu_table(picrust_phyloseq_cc)),
                              phyloseq::sample_data(picrust_phyloseq_cc)$Habitat,
                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)

# Create data frame and format
cc_ec_otu_table <- data.frame(phyloseq::otu_table(picrust_phyloseq_cc))
cc_ec_otu_table <- rownames_to_column(cc_ec_otu_table, var = "OTU")

write.csv(aldex2_ec_cc, "Picrust/Aldex2/aldex2_ec_cc.csv")

# Import results CSV and format
aldex2_ec_cc_result <- read.csv("Picrust/Aldex2/aldex2_ec_cc.csv")
colnames(aldex2_ec_cc_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

# Create data frame for taxonomy 
aldex_taxa_info_ec_cc <- data.frame(tax_table(picrust_phyloseq_cc))
aldex_taxa_info_ec_cc <- aldex_taxa_info_ec_cc %>%
  rownames_to_column(var = "OTU")

sample_tab_ec_cc <- data.frame(sample_data(picrust_phyloseq_cc))

# Filter aldex2 results by sig kw.ep and join the taxanomic information
sig_aldex2_ec_cc_result <- aldex2_ec_cc_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)

#sig_aldex2_ec_cc_result <- sig_aldex2_ec_cc_result %>% top_n(-20, kw.ep)

sig_aldex2_ec_cc_result <- left_join(sig_aldex2_ec_cc_result, aldex_taxa_info_ec_cc)

write.csv(sig_aldex2_ec_cc_result, "Picrust/Aldex2/sig_aldex2_ec_cc_result.csv")


# Made changes to taxonomy names if needed 
sig_aldex2_ec_cc_result <- read.csv("Picrust/Aldex2/sig_aldex2_ec_cc_result.csv")

sig_aldex2_ec_cc_result <- sig_aldex2_ec_cc_result[grepl("sulf", sig_aldex2_ec_cc_result$description, ignore.case = TRUE), ]


# Create clr objects by using OTU ids, original otu table, and all significant OTUs
sig_aldex2_ec_cc_result_count <- left_join(sig_aldex2_ec_cc_result, cc_ec_otu_table)
#sig_aldex2_ec_cc_result_count <- sig_aldex2_ec_cc_result_count[, -1]
write.csv(sig_aldex2_ec_cc_result_count, "Picrust/Aldex2/sig_aldex2_ec_cc_result_count_sulfur.csv")

# Dropped taxonomic rank and formatted columns 
clr_ec_cc <- sig_aldex2_ec_cc_result_count[, -(2:4)] 
rownames(clr_ec_cc) <- clr_ec_cc$OTU
clr_ec_cc <- clr_ec_cc[, -1]

# Adjusting zeros on the matrix and applying log transformation
shsk_ec_cc_czm <- cmultRepl(t(clr_ec_cc),  label=0, method="CZM")
shsk_ec_cc_czm_tv <- t(apply(shsk_ec_cc_czm, 1, function(x){log(x) - mean(log(x))}))
shsk_ec_cc_czm <- (apply(clr_ec_cc, 1, function(x){log(x+1) - mean(log(x+1))}))


# Combine the heatmap and the annotation
Z.Score.ec_cc <- scale(t(shsk_ec_cc_czm))

# Organize total abundance and taxa name
# Taxa name
Z.Score.ec_name_cc <- as.data.frame(Z.Score.ec_cc)
str(Z.Score.ec_name_cc)
Z.Score.ec_name_cc <- rownames_to_column(Z.Score.ec_name_cc, var = "OTU")
Z.Score.ec_name_cc <- left_join(Z.Score.ec_name_cc, aldex_taxa_info_ec_cc)
head(Z.Score.ec_name_cc)

# Taxa abundance 
ec_otu_table_total_cc <- as.data.frame(cc_ec_otu_table)
ec_otu_table_total_cc$Total <- rowSums(cc_ec_otu_table[, -1])
head(ec_otu_table_total_cc)
ec_otu_table_total_cc <- ec_otu_table_total_cc[, -(2:24)] 

Z.Score.ec_count_total_cc <- left_join(Z.Score.ec_name_cc, ec_otu_table_total_cc)
head(Z.Score.ec_count_total_cc)

# Defining color scheme for row annotations
abundance_col_fun_ec_cc = colorRamp2(c(0, 50000, 100000, 150000),
                                     c("#c7eae5",
                                       "#80cdc1",
                                       "#35978f",
                                       "#01665e"))

ha_right_ec_cc = rowAnnotation(
  Abundance = Z.Score.ec_count_total_cc$Total, border = FALSE, col = list(Abundance = abundance_col_fun_ec_cc))
row_labels_ec_cc = Z.Score.ec_count_total_cc$description


ha_ec_cc = HeatmapAnnotation(
  Habitat = as.vector(sample_tab_ec_cc$Habitat),
  col = list(
    Habitat = c("Bare Sediment" = "saddlebrown", "Sea Grass" = "#00A572")),
  annotation_legend_param = list(
    Habitat = list(
      title = "Habitat",
      at = c("Bare Sediment", "Sea Grass"),
      labels = c ("BS", "SG")
    )
  ))

# Plot heatmap at the EC level
hm_ec_cc <- Heatmap(Z.Score.ec_cc, name = "Z-score, CLR", col = col_matrix,
                    column_title  = "CC Sulfur EC", 
                    column_title_gp = gpar(fontface = "bold", fontsize = 14),
                    #column_split = as.vector(as.vector(sample_tab_ec_cc$Habitat)),
                    #column_order = order(as.numeric(gsub("column", "", colnames(Z.Score.gen)))),
                    border = TRUE,
                    top_annotation = ha_ec_cc,
                    right_annotation = ha_right_ec_cc,
                    row_title = "EC",
                    row_labels = row_labels_ec_cc,
                    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                    row_names_gp = gpar(fontsize = 6),
                    column_names_gp = gpar(fontsize = 6),
                    #row_order = order(row_labels_ec_cc),
                    rect_gp = gpar(col = "white", lwd = 1),
                    show_column_names = FALSE,
                    show_heatmap_legend = TRUE)

hm_ec_cc


####################################### Run ALEDx2 for Westside Park only ###################################
picrust_phyloseq_wp <- subset_samples(picrust_phyloseq, Site == "Westside Park")

d1 <- psmelt(picrust_phyloseq_wp)

aldex2_ec_wp <- ALDEx2::aldex(data.frame(phyloseq::otu_table(picrust_phyloseq_wp)),
                              phyloseq::sample_data(picrust_phyloseq_wp)$Habitat,
                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)

# Create data frame and format
wp_ec_otu_table <- data.frame(phyloseq::otu_table(picrust_phyloseq_wp))
wp_ec_otu_table <- rownames_to_column(wp_ec_otu_table, var = "OTU")

write.csv(aldex2_ec_wp, "Picrust/Aldex2/aldex2_ec_wp.csv")

# Import results CSV and format
aldex2_ec_wp_result <- read.csv("Picrust/Aldex2/aldex2_ec_wp.csv")
colnames(aldex2_ec_wp_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

# Create data frame for taxonomy 
aldex_taxa_info_ec_wp <- data.frame(tax_table(picrust_phyloseq_wp))
aldex_taxa_info_ec_wp <- aldex_taxa_info_ec_wp %>%
  rownames_to_column(var = "OTU")

sample_tab_ec_wp <- data.frame(sample_data(picrust_phyloseq_wp))

# Filter aldex2 results by sig kw.ep and join the taxanomic information
sig_aldex2_ec_wp_result <- aldex2_ec_wp_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)

#sig_aldex2_ec_wp_result <- sig_aldex2_ec_wp_result %>% top_n(-20, kw.ep)

sig_aldex2_ec_wp_result <- left_join(sig_aldex2_ec_wp_result, aldex_taxa_info_ec_wp)

write.csv(sig_aldex2_ec_wp_result, "Picrust/Aldex2/sig_aldex2_ec_wp_result.csv")


# Made changes to taxonomy names if needed 
sig_aldex2_ec_wp_result <- read.csv("Picrust/Aldex2/sig_aldex2_ec_wp_result.csv")

sig_aldex2_ec_wp_result <- sig_aldex2_ec_wp_result[grepl("sulf", sig_aldex2_ec_wp_result$description, ignore.case = TRUE), ]


# Create clr objects by using OTU ids, original otu table, and all significant OTUs
sig_aldex2_ec_wp_result_count <- left_join(sig_aldex2_ec_wp_result, wp_ec_otu_table)
sig_aldex2_ec_wp_result_count <- sig_aldex2_ec_wp_result_count[, -1]
write.csv(sig_aldex2_ec_wp_result_count, "Picrust/Aldex2/sig_aldex2_ec_wp_result_count_sulfur.csv")

# Dropped taxonomic rank and formatted columns 
clr_ec_wp <- sig_aldex2_ec_wp_result_count[, -(2:4)] 
rownames(clr_ec_wp) <- clr_ec_wp$OTU
clr_ec_wp <- clr_ec_wp[, -1]

# Adjusting zeros on the matrix and applying log transformation
shsk_ec_wp_czm <- cmultRepl(t(clr_ec_wp),  label=0, method="CZM")
shsk_ec_wp_czm_tv <- t(apply(shsk_ec_wp_czm, 1, function(x){log(x) - mean(log(x))}))
shsk_ec_wp_czm <- (apply(clr_ec_wp, 1, function(x){log(x+1) - mean(log(x+1))}))


# Combine the heatmap and the annotation
Z.Score.ec_wp <- scale(t(shsk_ec_wp_czm))

# Organize total abundance and taxa name
# Taxa name
Z.Score.ec_name_wp <- as.data.frame(Z.Score.ec_wp)
str(Z.Score.ec_name_wp)
Z.Score.ec_name_wp <- rownames_to_column(Z.Score.ec_name_wp, var = "OTU")
Z.Score.ec_name_wp <- left_join(Z.Score.ec_name_wp, aldex_taxa_info_ec_wp)
head(Z.Score.ec_name_wp)

# Taxa abundance 
ec_otu_table_total_wp <- as.data.frame(wp_ec_otu_table)
ec_otu_table_total_wp$Total <- rowSums(wp_ec_otu_table[, -1])
head(ec_otu_table_total_wp)
ec_otu_table_total_wp <- ec_otu_table_total_wp[, -(2:37)] 

Z.Score.ec_count_total_wp <- left_join(Z.Score.ec_name_wp, ec_otu_table_total_wp)
head(Z.Score.ec_count_total_wp)

# Defining color scheme for row annotations
abundance_col_fun_ec_wp = colorRamp2(c(0, 200000, 500000, 1000000),
                                     c("#c7eae5",
                                       "#80cdc1",
                                       "#35978f",
                                       "#01665e"))

ha_right_ec_wp = rowAnnotation(
  Abundance = Z.Score.ec_count_total_wp$Total, border = FALSE, col = list(Abundance = abundance_col_fun_ec_wp))
row_labels_ec_wp = Z.Score.ec_count_total_wp$description


ha_ec_wp = HeatmapAnnotation(
  Habitat = as.vector(sample_tab_ec_wp$Habitat),
  col = list(
    Habitat = c("Bare Sediment" = "saddlebrown", "Sea Grass" = "#00A572")),
  annotation_legend_param = list(
    Habitat = list(
      title = "Habitat",
      at = c("Bare Sediment", "Sea Grass"),
      labels = c ("BS", "SG")
    )
  ))

# Plot heatmap at the EC level
hm_ec_wp <- Heatmap(Z.Score.ec_wp, name = "Z-score, CLR", col = col_matrix,
                    column_title  = "WP Sulfur EC", 
                    column_title_gp = gpar(fontface = "bold", fontsize = 14),
                    #column_split = as.vector(as.vector(sample_tab_ec_wp$Habitat)),
                    #column_order = order(as.numeric(gsub("column", "", colnames(Z.Score.gen)))),
                    border = TRUE,
                    top_annotation = ha_ec_wp,
                    right_annotation = ha_right_ec_wp,
                    row_title = "EC",
                    row_labels = row_labels_ec_wp,
                    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                    row_names_gp = gpar(fontsize = 6),
                    column_names_gp = gpar(fontsize = 6),
                    #row_order = order(row_labels_ec_wp),
                    rect_gp = gpar(col = "white", lwd = 1),
                    show_column_names = FALSE,
                    show_heatmap_legend = TRUE)

hm_ec_wp

####################################### Run ALEDx2 for Mason's Marina only ###################################
picrust_phyloseq_mm <- subset_samples(picrust_phyloseq, Site == "Mason's Marina")

d1 <- psmelt(picrust_phyloseq_mm)

aldex2_ec_mm <- ALDEx2::aldex(data.frame(phyloseq::otu_table(picrust_phyloseq_mm)),
                              phyloseq::sample_data(picrust_phyloseq_mm)$Habitat,
                              mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)

# Create data frame and format
mm_ec_otu_table <- data.frame(phyloseq::otu_table(picrust_phyloseq_mm))
mm_ec_otu_table <- rownames_to_column(mm_ec_otu_table, var = "OTU")

write.csv(aldex2_ec_mm, "Picrust/Aldex2/aldex2_ec_mm.csv")

# Import results CSV and format
aldex2_ec_mm_result <- read.csv("Picrust/Aldex2/aldex2_ec_mm.csv")
colnames(aldex2_ec_mm_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

# Create data frame for taxonomy 
aldex_taxa_info_ec_mm <- data.frame(tax_table(picrust_phyloseq_mm))
aldex_taxa_info_ec_mm <- aldex_taxa_info_ec_mm %>%
  rownames_to_column(var = "OTU")

sample_tab_ec_mm <- data.frame(sample_data(picrust_phyloseq_mm))

# Filter aldex2 results by sig kw.ep and join the taxanomic information
sig_aldex2_ec_mm_result <- aldex2_ec_mm_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)

#sig_aldex2_ec_mm_result <- sig_aldex2_ec_mm_result %>% top_n(-20, kw.ep)

sig_aldex2_ec_mm_result <- left_join(sig_aldex2_ec_mm_result, aldex_taxa_info_ec_mm)

write.csv(sig_aldex2_ec_mm_result, "Picrust/Aldex2/sig_aldex2_ec_mm_result.csv")


# Made changes to taxonomy names if needed 
sig_aldex2_ec_mm_result <- read.csv("Picrust/Aldex2/sig_aldex2_ec_mm_result.csv")

sig_aldex2_ec_mm_result <- sig_aldex2_ec_mm_result[grepl("sulf", sig_aldex2_ec_mm_result$description, ignore.case = TRUE), ]

# Create clr objects by using OTU ids, original otu table, and all significant OTUs
sig_aldex2_ec_mm_result_count <- left_join(sig_aldex2_ec_mm_result, mm_ec_otu_table)
sig_aldex2_ec_mm_result_count <- sig_aldex2_ec_mm_result_count[, -1]
write.csv(sig_aldex2_ec_mm_result_count, "Picrust/Aldex2/sig_aldex2_ec_mm_result_count.csv")

# Dropped taxonomic rank and formatted columns 
clr_ec_mm <- sig_aldex2_ec_mm_result_count[, -(2:4)] 
rownames(clr_ec_mm) <- clr_ec_mm$OTU
clr_ec_mm <- clr_ec_mm[, -1]

# Adjusting zeros on the matrix and applying log transformation
shsk_ec_mm_czm <- cmultRepl(t(clr_ec_mm),  label=0, method="CZM")
shsk_ec_mm_czm_tv <- t(apply(shsk_ec_mm_czm, 1, function(x){log(x) - mean(log(x))}))
shsk_ec_mm_czm <- (apply(clr_ec_mm, 1, function(x){log(x+1) - mean(log(x+1))}))


# Combine the heatmap and the annotation
Z.Score.ec_mm <- scale(t(shsk_ec_mm_czm))

# Organize total abundance and taxa name
# Taxa name
Z.Score.ec_name_mm <- as.data.frame(Z.Score.ec_mm)
str(Z.Score.ec_name_mm)
Z.Score.ec_name_mm <- rownames_to_column(Z.Score.ec_name_mm, var = "OTU")
Z.Score.ec_name_mm <- left_join(Z.Score.ec_name_mm, aldex_taxa_info_ec_mm)
head(Z.Score.ec_name_mm)

# Taxa abundance 
ec_otu_table_total_mm <- as.data.frame(mm_ec_otu_table)
ec_otu_table_total_mm$Total <- rowSums(mm_ec_otu_table[, -1])
head(ec_otu_table_total_mm)
ec_otu_table_total_mm <- ec_otu_table_total_mm[, -(2:24)] 

Z.Score.ec_count_total_mm <- left_join(Z.Score.ec_name_mm, ec_otu_table_total_mm)
head(Z.Score.ec_count_total_mm)

# Defining color scheme for row annotations
abundance_col_fun_ec_mm = colorRamp2(c(0, 10000, 20000, 30000),
                                     c("#c7eae5",
                                       "#80cdc1",
                                       "#35978f",
                                       "#01665e"))

ha_right_ec_mm = rowAnnotation(
  Abundance = Z.Score.ec_count_total_mm$Total, border = FALSE, col = list(Abundance = abundance_col_fun_ec_mm))
row_labels_ec_mm = Z.Score.ec_count_total_mm$description


ha_ec_mm = HeatmapAnnotation(
  Habitat = as.vector(sample_tab_ec_mm$Habitat),
  col = list(
    Habitat = c("Bare Sediment" = "saddlebrown", "Sea Grass" = "#00A572")),
  annotation_legend_param = list(
    Habitat = list(
      title = "Habitat",
      at = c("Bare Sediment", "Sea Grass"),
      labels = c ("BS", "SG")
    )
  ))

# Plot heatmap at the EC level
hm_ec_mm <- Heatmap(Z.Score.ec_mm, name = "Z-score, CLR", col = col_matrix,
                    column_title  = "MM Sulfur EC", 
                    column_title_gp = gpar(fontface = "bold", fontsize = 14),
                    #column_split = as.vector(as.vector(sample_tab_ec_mm$Habitat)),
                    #column_order = order(as.numeric(gsub("column", "", colnames(Z.Score.gen)))),
                    border = TRUE,
                    top_annotation = ha_ec_mm,
                    right_annotation = ha_right_ec_mm,
                    row_title = "EC",
                    row_labels = row_labels_ec_mm,
                    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                    row_names_gp = gpar(fontsize = 6),
                    column_names_gp = gpar(fontsize = 6),
                    #row_order = order(row_labels_ec_mm),
                    rect_gp = gpar(col = "white", lwd = 1),
                    show_column_names = FALSE,
                    show_heatmap_legend = TRUE)

hm_ec_mm

####################################### Run ALEDx2 for Bare Sediment only ###################################
picrust_phyloseq_bare <- subset_samples(picrust_phyloseq, Habitat == "Bare Sediment")

d1 <- psmelt(picrust_phyloseq_bare)

aldex2_ec_bare <- ALDEx2::aldex(data.frame(phyloseq::otu_table(picrust_phyloseq_bare)),
                                phyloseq::sample_data(picrust_phyloseq_bare)$Site,
                                mc.samples=128, test="kw", effect=TRUE,include.sample.subareary=FALSE, denom="all", verbose=FALSE)

# Create data frame and format
bare_ec_otu_table <- data.frame(phyloseq::otu_table(picrust_phyloseq_bare))
bare_ec_otu_table <- rownames_to_column(bare_ec_otu_table, var = "OTU")

write.csv(aldex2_ec_bare, "Picrust/Aldex2/aldex2_ec_bare.csv")

# Import results CSV and format
aldex2_ec_bare_result <- read.csv("Picrust/Aldex2/aldex2_ec_bare.csv")
colnames(aldex2_ec_bare_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

# Create data frame for taxonomy 
aldex_taxa_info_ec_bare <- data.frame(tax_table(picrust_phyloseq_bare))
aldex_taxa_info_ec_bare <- aldex_taxa_info_ec_bare %>%
  rownames_to_column(var = "OTU")

sample_tab_ec_bare <- data.frame(sample_data(picrust_phyloseq_bare))

# Filter aldex2 results by sig kw.ep and join the taxanomic information
sig_aldex2_ec_bare_result <- aldex2_ec_bare_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)

sig_aldex2_ec_bare_result <- left_join(sig_aldex2_ec_bare_result, aldex_taxa_info_ec_bare)

write.csv(sig_aldex2_ec_bare_result, "Picrust/Aldex2/sig_aldex2_ec_bare_result.csv")


# Made changes to taxonomy names if needed 
sig_aldex2_ec_bare_result <- read.csv("Picrust/Aldex2/sig_aldex2_ec_bare_result.csv")

sig_aldex2_ec_bare_result <- sig_aldex2_ec_bare_result[grepl("sulf", sig_aldex2_ec_bare_result$description, ignore.case = TRUE), ]

# Create clr objects by using OTU ids, original otu table, and all significant OTUs
sig_aldex2_ec_bare_result_count <- left_join(sig_aldex2_ec_bare_result, bare_ec_otu_table)
sig_aldex2_ec_bare_result_count <- sig_aldex2_ec_bare_result_count[, -1]
write.csv(sig_aldex2_ec_bare_result_count, "Picrust/Aldex2/sig_aldex2_ec_bare_result_count.csv")

# Dropped taxonomic rank and formatted columns 
clr_ec_bare <- sig_aldex2_ec_bare_result_count[, -(2:4)] 
rownames(clr_ec_bare) <- clr_ec_bare$OTU
clr_ec_bare <- clr_ec_bare[, -1]

# Adjusting zeros on the matrix and applying log transformation
shsk_ec_bare_czm <- cmultRepl(t(clr_ec_bare),  label=0, method="CZM")
shsk_ec_bare_czm_tv <- t(apply(shsk_ec_bare_czm, 1, function(x){log(x) - mean(log(x))}))
shsk_ec_bare_czm <- (apply(clr_ec_bare, 1, function(x){log(x+1) - mean(log(x+1))}))


# Combine the heatmap and the annotation
Z.Score.ec_bare <- scale(t(shsk_ec_bare_czm))

# Organize total abundance and taxa name
# Taxa name
Z.Score.ec_name_bare <- as.data.frame(Z.Score.ec_bare)
str(Z.Score.ec_name_bare)
Z.Score.ec_name_bare <- rownames_to_column(Z.Score.ec_name_bare, var = "OTU")
Z.Score.ec_name_bare <- left_join(Z.Score.ec_name_bare, aldex_taxa_info_ec_bare)
head(Z.Score.ec_name_bare)

# Taxa abundance 
ec_otu_table_total_bare <- as.data.frame(bare_ec_otu_table)
ec_otu_table_total_bare$Total <- rowSums(bare_ec_otu_table[, -1])
head(ec_otu_table_total_bare)
ec_otu_table_total_bare <- ec_otu_table_total_bare[, -(2:49)] 

Z.Score.ec_count_total_bare <- left_join(Z.Score.ec_name_bare, ec_otu_table_total_bare)
head(Z.Score.ec_count_total_bare)

# Defining color scheme for row annotations
abundance_col_fun_ec_bare = colorRamp2(c(0, 200000, 500000, 1000000),
                                       c("#c7eae5",
                                         "#80cdc1",
                                         "#35978f",
                                         "#01665e"))

ha_right_ec_bare = rowAnnotation(
  Abundance = Z.Score.ec_count_total_bare$Total, border = FALSE, col = list(Abundance = abundance_col_fun_ec_bare))
row_labels_ec_bare = Z.Score.ec_count_total_bare$description

ha_ec_bare = HeatmapAnnotation(
  Site = as.vector(sample_tab_ec_bare$Site),
  col = list(
    Site = c("Mason's Marina" = "#E43F6F", "Westside Park" = "#2E5266",
             "Campbell Cove" = "#C8AB83")
  ),
  Site = list(
    title = "Site",
    at = c("Campbell Cove", "Westside Park", "Mason's Marina"),
    labels = c ("CC", "WP", "MM")
  )
)

# Plot heatmap at the EC level
hm_ec_bare <- Heatmap(Z.Score.ec_bare, name = "Z-score, CLR", col = col_matrix,
                      column_title  = "Bare Sulf", 
                      column_title_gp = gpar(fontface = "bold", fontsize = 14),
                      #column_split = as.vector(as.vector(sample_tab_ec_bare$Site)),
                      #column_order = order(as.numeric(gsub("column", "", colnames(Z.Score.gen)))),
                      border = TRUE,
                      top_annotation = ha_ec_bare,
                      right_annotation = ha_right_ec_bare,
                      row_title = "EC",
                      row_labels = row_labels_ec_bare,
                      row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                      row_names_gp = gpar(fontsize = 6),
                      column_names_gp = gpar(fontsize = 6),
                      #row_order = order(row_labels_ec_bare),
                      rect_gp = gpar(col = "white", lwd = 1),
                      show_column_names = FALSE,
                      show_heatmap_legend = TRUE)

hm_ec_bare


####################################### Run ALEDx2 for Sea Grass only ###################################

picrust_phyloseq_sg <- subset_samples(picrust_phyloseq, Habitat == "Sea Grass")

d1 <- psmelt(picrust_phyloseq_sg)

aldex2_ec_sg <- ALDEx2::aldex(data.frame(phyloseq::otu_table(picrust_phyloseq_sg)),
                              phyloseq::sample_data(picrust_phyloseq_sg)$Site,
                              mc.samples=128, test="kw", effect=TRUE,include.sample.susgary=FALSE, denom="all", verbose=FALSE)

# Create data frame and format
sg_ec_otu_table <- data.frame(phyloseq::otu_table(picrust_phyloseq_sg))
sg_ec_otu_table <- rownames_to_column(sg_ec_otu_table, var = "OTU")

write.csv(aldex2_ec_sg, "Picrust/Aldex2/aldex2_ec_sg.csv")

# Import results CSV and format
aldex2_ec_sg_result <- read.csv("Picrust/Aldex2/aldex2_ec_sg.csv")
colnames(aldex2_ec_sg_result) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

# Create data frame for taxonomy 
aldex_taxa_info_ec_sg <- data.frame(tax_table(picrust_phyloseq_sg))
aldex_taxa_info_ec_sg <- aldex_taxa_info_ec_sg %>%
  rownames_to_column(var = "OTU")

sample_tab_ec_sg <- data.frame(sample_data(picrust_phyloseq_sg))

# Filter aldex2 results by sig kw.ep and join the taxanomic information
sig_aldex2_ec_sg_result <- aldex2_ec_sg_result %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)

sig_aldex2_ec_sg_result <- left_join(sig_aldex2_ec_sg_result, aldex_taxa_info_ec_sg)

write.csv(sig_aldex2_ec_sg_result, "Picrust/Aldex2/sig_aldex2_ec_sg_result.csv")


# Made changes to taxonomy names if needed 
sig_aldex2_ec_sg_result <- read.csv("Picrust/Aldex2/sig_aldex2_ec_sg_result.csv")

sig_aldex2_ec_sg_result <- sig_aldex2_ec_sg_result[grepl("sulf", sig_aldex2_ec_sg_result$description, ignore.case = TRUE), ]

# Create clr objects by using OTU ids, original otu table, and all significant OTUs
sig_aldex2_ec_sg_result_count <- left_join(sig_aldex2_ec_sg_result, sg_ec_otu_table)
sig_aldex2_ec_sg_result_count <- sig_aldex2_ec_sg_result_count[, -1]
write.csv(sig_aldex2_ec_sg_result_count, "Picrust/Aldex2/sig_aldex2_ec_sg_result_count.csv")

# Dropped taxonomic rank and formatted columns 
clr_ec_sg <- sig_aldex2_ec_sg_result_count[, -(2:4)] 
rownames(clr_ec_sg) <- clr_ec_sg$OTU
clr_ec_sg <- clr_ec_sg[, -1]

# Adjusting zeros on the matrix and applying log transformation
shsk_ec_sg_czm <- cmultRepl(t(clr_ec_sg),  label=0, method="CZM")
shsk_ec_sg_czm_tv <- t(apply(shsk_ec_sg_czm, 1, function(x){log(x) - mean(log(x))}))
shsk_ec_sg_czm <- (apply(clr_ec_sg, 1, function(x){log(x+1) - mean(log(x+1))}))


# Combine the heatmap and the annotation
Z.Score.ec_sg <- scale(t(shsk_ec_sg_czm))

# Organize total abundance and taxa name
# Taxa name
Z.Score.ec_name_sg <- as.data.frame(Z.Score.ec_sg)
str(Z.Score.ec_name_sg)
Z.Score.ec_name_sg <- rownames_to_column(Z.Score.ec_name_sg, var = "OTU")
Z.Score.ec_name_sg <- left_join(Z.Score.ec_name_sg, aldex_taxa_info_ec_sg)
head(Z.Score.ec_name_sg)

# Taxa abundance 
ec_otu_table_total_sg <- as.data.frame(sg_ec_otu_table)
ec_otu_table_total_sg$Total <- rowSums(sg_ec_otu_table[, -1])
head(ec_otu_table_total_sg)
ec_otu_table_total_sg <- ec_otu_table_total_sg[, -(2:37)] 

Z.Score.ec_count_total_sg <- left_join(Z.Score.ec_name_sg, ec_otu_table_total_sg)
head(Z.Score.ec_count_total_sg)

# Defining color scheme for row annotations
abundance_col_fun_ec_sg = colorRamp2(c(0, 50000, 150000, 500000),
                                     c("#c7eae5",
                                       "#80cdc1",
                                       "#35978f",
                                       "#01665e"))

ha_right_ec_sg = rowAnnotation(
  Abundance = Z.Score.ec_count_total_sg$Total, border = FALSE, col = list(Abundance = abundance_col_fun_ec_sg))
row_labels_ec_sg = Z.Score.ec_count_total_sg$description

ha_ec_sg = HeatmapAnnotation(
  Site = as.vector(sample_tab_ec_sg$Site),
  col = list(
    Site = c("Mason's Marina" = "#E43F6F", "Westside Park" = "#2E5266",
             "Campbell Cove" = "#C8AB83")
  ),
  Site = list(
    title = "Site",
    at = c("Campbell Cove", "Westside Park", "Mason's Marina"),
    labels = c ("CC", "WP", "MM")
  )
)

# Plot heatmap at the EC level
hm_ec_sg <- Heatmap(Z.Score.ec_sg, name = "Z-score, CLR", col = col_matrix,
                    column_title  = "Sea Grass Sulf", 
                    column_title_gp = gpar(fontface = "bold", fontsize = 14),
                    column_split = as.vector(as.vector(sample_tab_ec_sg$Site)),
                    #column_order = order(as.numeric(gsub("column", "", colnames(Z.Score.gen)))),
                    border = TRUE,
                    top_annotation = ha_ec_sg,
                    right_annotation = ha_right_ec_sg,
                    row_title = "EC",
                    row_labels = row_labels_ec_sg,
                    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                    row_names_gp = gpar(fontsize = 6),
                    column_names_gp = gpar(fontsize = 6),
                    #row_order = order(row_labels_ec_sg),
                    rect_gp = gpar(col = "white", lwd = 1),
                    show_column_names = FALSE,
                    show_heatmap_legend = TRUE)

hm_ec_sg

####################################### Run ALEDx2 for habitat ###################################
aldex2_ec_habitat <- ALDEx2::aldex(data.frame(phyloseq::otu_table(picrust_phyloseq)),
                           phyloseq::sample_data(picrust_phyloseq)$Habitat,
                           mc.samples=128, test="kw", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
# Create data frame and format
ec_otu_table <- data.frame(phyloseq::otu_table(picrust_phyloseq))
ec_otu_table <- rownames_to_column(ec_otu_table, var = "OTU")

write.csv(aldex2_ec_habitat, "Picrust/Aldex2/aldex2_ec_habitat.csv")

# Import results CSV and format
aldex2_ec_result_habitat <- read.csv("Picrust/Aldex2/aldex2_ec_habitat.csv")
colnames(aldex2_ec_result_habitat) <- c("OTU", "kw.ep", "kw.eBH", "glm.ep", "glm.eBH")

# Create data frame for taxonomy 
aldex_taxa_info_ec <- data.frame(tax_table(picrust_phyloseq))
aldex_taxa_info_ec <- aldex_taxa_info_ec %>%
  rownames_to_column(var = "OTU")

sample_tab_ec <- data.frame(sample_data(picrust_phyloseq))

# Filter aldex2 results by sig kw.ep and join the taxanomic information
sig_aldex2_ec_result_habitat <- aldex2_ec_result_habitat %>%
  filter(kw.ep < 0.05) %>%
  arrange(kw.ep, kw.eBH) %>%
  dplyr::select(OTU, kw.ep, kw.eBH)

sig_aldex2_ec_result_habitat <- left_join(sig_aldex2_ec_result_habitat, aldex_taxa_info_ec)


write.csv(sig_aldex2_ec_result_habitat, "Picrust/Aldex2/sig_aldex2_ec_result_habitat.csv")


# Made changes to taxonomy names if needed 
sig_aldex2_ec_result_habitat <- read.csv("Picrust/Aldex2/sig_aldex2_ec_result_habitat.csv")

sig_aldex2_ec_result_habitat <- sig_aldex2_ec_result_habitat[grepl("sulf", sig_aldex2_ec_result_habitat$description, ignore.case = TRUE), ]


# Create clr objects by using OTU ids, original otu table, and all significant OTUs
sig_aldex2_ec_result_count_habitat <- left_join(sig_aldex2_ec_result_habitat, ec_otu_table)
sig_aldex2_ec_result_count_habitat <- sig_aldex2_ec_result_count_habitat[, -1]
write.csv(sig_aldex2_ec_result_count_habitat, "Picrust/Aldex2/sig_aldex2_ec_result_count_habitat.csv")

# Dropped taxonomic rank and formatted columns 
clr_ec_hab <- sig_aldex2_ec_result_count_habitat[, -(2:4)] 
rownames(clr_ec_hab) <- clr_ec_hab$OTU
clr_ec_hab <- clr_ec_hab[, -1]

# Adjusting zeros on the matrix and applying log transformation
shsk_ec_czm_hab <- cmultRepl(t(clr_ec_hab),  label=0, method="CZM")
shsk_ec_czm_tv_hab <- t(apply(shsk_ec_czm_hab, 1, function(x){log(x) - mean(log(x))}))
shsk_ec_czm_hab <- (apply(clr_ec_hab, 1, function(x){log(x+1) - mean(log(x+1))}))


# Combine the heatmap and the annotation
Z.Score.ec_hab <- scale(t(shsk_ec_czm_hab))

# Organize total abundance and taxa name
# Taxa name
Z.Score.ec_name_hab <- as.data.frame(Z.Score.ec_hab)
str(Z.Score.ec_name_hab)
Z.Score.ec_name_hab <- rownames_to_column(Z.Score.ec_name_hab, var = "OTU")
Z.Score.ec_name_hab <- left_join(Z.Score.ec_name_hab, aldex_taxa_info_ec)
head(Z.Score.ec_name_hab)

# Taxa abundance 
ec_otu_table_total_hab <- as.data.frame(ec_otu_table)
ec_otu_table_total_hab$Total <- rowSums(ec_otu_table_total_hab[, -1])
head(ec_otu_table_total_hab)
ec_otu_table_total_hab <- ec_otu_table_total_hab[, -(2:85)] 

Z.Score.ec_count_total_hab <- left_join(Z.Score.ec_name_hab, ec_otu_table_total_hab)
head(Z.Score.ec_count_total)

# Defining color scheme for row annotations
abundance_col_fun_ec_hab = colorRamp2(c(0, 200000, 500000, 1000000),
                                  c("#c7eae5",
                                    "#80cdc1",
                                    "#35978f",
                                    "#01665e"))

ha_right_ec_hab = rowAnnotation(
  Abundance = Z.Score.ec_count_total_hab$Total, border = FALSE, col = list(Abundance = abundance_col_fun_ec_hab))
row_labels_ec_hab = Z.Score.ec_count_total_hab$description


ha_ec = HeatmapAnnotation(
  Habitat = as.vector(sample_tab_ec$Habitat),
  Site = as.vector(sample_tab_ec$Site),
  col = list(
    Habitat = c("Bare Sediment" = "saddlebrown", "Sea Grass" = "#00A572"),
    Site = c("Mason's Marina" = "#E43F6F", "Westside Park" = "#2E5266",
             "Campbell Cove" = "#C8AB83")
  ),
  annotation_legend_param = list(
    Habitat = list(
      title = "Habitat",
      at = c("Bare Sediment", "Sea Grass"),
      labels = c ("BS", "SG")
    ),
    Site = list(
      title = "Site",
      at = c("Campbell Cove", "Westside Park", "Mason's Marina"),
      labels = c ("CC", "WP", "MM")
    )
  ))

# Plot heatmap at the EC level
hm_ec_hab <- Heatmap(Z.Score.ec_hab, name = "Z-score, CLR", col = col_matrix,
                 column_title  = "Sites and Habitats", 
                 column_title_gp = gpar(fontface = "bold", fontsize = 14),
                 column_split = as.vector(as.vector(sample_tab_ec$Site)),
                 column_order = order(as.numeric(gsub("column", "", colnames(Z.Score.ec_hab)))),
                 border = TRUE,
                 top_annotation = ha_ec,
                 right_annotation = ha_right_ec_hab,
                 row_title = "EC",
                 row_labels = row_labels_ec_hab,
                 row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 row_names_gp = gpar(fontsize = 6),
                 column_names_gp = gpar(fontsize = 6),
                 #row_order = order(row_labels_ec),
                 rect_gp = gpar(col = "white", lwd = 1),
                 show_column_names = FALSE,
                 show_heatmap_legend = TRUE)


hm_ec_hab
