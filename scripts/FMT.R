pacman::p_load(RColorBrewer, phyloseq, tidyverse, patchwork, here, ggplot2)

#Loading dataset ----
fmt_metagenomes <- import_biom("data/FMT paper merged biom.biom2")
view(fmt_metagenomes@sam_data@.Data)

# inspect objects ----
class(fmt_metagenomes) # check the class of the imported object
sample_names(fmt_metagenomes) #check the samples names in the phyloseq object
view(fmt_metagenomes@sam_data@.Data) # view the sample data
view(fmt_metagenomes@tax_table) # view the taxonomy table
view(fmt_metagenomes@otu_table) # view the OTU table 
otu_table(fmt_metagenomes) # Get OTU/sample counts as a matrix

# Show the most abundant taxa in each sample (at species level)
top_species <- apply(otu_table(fmt_metagenomes), 2, function(x) {
  top <- sort(x, decreasing = TRUE)[1:5]  # top 5 taxa
  return(names(top))
})

# This prints a matrix: each column = sample, each row = top taxa
print(top_species)

# Preview the sample names ----
sample_names(fmt_metagenomes)

# Update sample names
sample_names(fmt_metagenomes) <- c("Donor4", "Donor3", "Donor2", "Donor1", "Non_responder4", "Non_responder3", "Non_responder2", "Non_responder1", "Responder4", "Responder3", "Responder2", "Responder1")

#Remove unnecessary characters in .data matrix ----
fmt_metagenomes@tax_table@.Data <- 
  substring(fmt_metagenomes@tax_table@.Data, 4)

#Rename the column names in the object
colnames(fmt_metagenomes@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Exploring the taxa present 
unique(fmt_metagenomes@tax_table@.Data[,"Kingdom"])
sum(fmt_metagenomes@tax_table@.Data[,"Phylum"] == "Proteobacteria")

# Subset and work with only the Bacteria kingdom
fmt_metagenomes@tax_table <- subset_taxa(fmt_metagenomes@tax_table, Kingdom == "Bacteria")
summary(fmt_metagenomes@tax_table@.Data == "")

#Diversity calculations ----
#alpha diversity calculation 
diversity_values <- 
  estimate_richness(fmt_metagenomes,
                    measures = c("Shannon", "Observed", "Chao1"))
diversity_values

# Visualizing the diversity plot
plot_richness(physeq = fmt_metagenomes,
              title = "Alpha diversity indexes for the twelve metagenome samples",
              measures = c("Shannon", "Observed", "Chao1"),
              sortby = "Shannon")

#beta diversity calculation
fmt_metagenomes_percentages <- 
  transform_sample_counts(fmt_metagenomes,function(x) x*100/sum(x))
head(fmt_metagenomes_percentages@otu_table@.Data)

#compare the abundances given by the percentage
meta_ord <- ordinate(physeq = fmt_metagenomes_percentages, method = "NMDS", distance = "bray")
plot_ordination(physeq = fmt_metagenomes_percentages, ordination = meta_ord)

# We can align the points to the metadata
metadata_fmt <- data.frame(Sample=c("Donor4", "Donor3", "Donor2", "Donor1", "Non_responder4", "Non_responder3", "Non_responder2", "Non_responder1", "Responder4", "Responder3", "Responder2", "Responder1"), Treatment=c("Stool donor", "Stool donor", "Stool donor", "Stool donor", "Recipient non-responder before FMT", "Recipient non-responder before FMT", "Recipient non-responder before FMT", "Recipient non-responder before FMT", "Recipient responder before FMT", "Recipient responder before FMT", "Recipient responder before FMT","Recipient responder before FMT")) # Making dataframe with metadata  

rownames(metadata_fmt) <- metadata_fmt$Sample # Using sample names as row names  
fmt_metagenomes_percentages@sam_data <- sample_data(metadata_fmt) # Adding metadata to sam_data table of phyloseq object percentages  
meta_ord <- ordinate(physeq = fmt_metagenomes_percentages, method = "NMDS", distance = "bray") # Calculating beta diversity    
plot_ordination(physeq = fmt_metagenomes_percentages, ordination = meta_ord, color = "Treatment", title = "Beta diversity indexes for the twelve metagenome samples") + # Plotting beta diversity.  
  geom_text(mapping = aes(label = Sample), size = 3, vjust = 1.5) 

# Taxonomy classification ----
# to make a more personalised visualization of the taxa in our data, we can group all the OTUs that have the same taxonomy at a certain taxonomic rank using the tax_glom function 
percentages_glom <- tax_glom(fmt_metagenomes_percentages, taxrank = 'Phylum')
View(percentages_glom@tax_table@.Data)

# The psmelt() function can also be used to melt phyloseq objects into a data.frame to manipulate them with packages like ggplot2 and vegan.
percentages_df <- psmelt(percentages_glom)
str(percentages_df)

# Create another dataframe with the original data, then compare relative and absolute abundances
absolute_glom <- tax_glom(physeq = fmt_metagenomes, taxrank = "Phylum")
absolute_df <- psmelt(absolute_glom)
str(absolute_df)

# Create a color palette with colorRamppalette that allows for the personalization of the plot. With colorRampPalette, we will choose eight colors from the Dark2 palette and make a “ramp” with it; that is, convert those eight colors to the number of colors needed to have one for each phylum in our data frame. We need to have our Phylum column in the factor structure for this.
absolute_df$Phylum <- as.factor(absolute_df$Phylum)
phylum_colors_abs<- colorRampPalette(brewer.pal(8,"Dark2")) (length(levels(absolute_df$Phylum)))

# Considering that there are too many taxa to adequately distinguish the color of each one, and less of the ones that hold the most incredible abundance. It is best to change the identification of the OTUs whose relative abundance is less than 0.2%:
percentages_df$Phylum <- as.character(percentages_df$Phylum) # Return the Phylum column to be of type character
percentages_df$Phylum[percentages_df$Abundance < 0.5] <- "Phyla < 0.5% abund."
unique(percentages_df$Phylum)

# Create the figure for the absolute data with absolute abundances (i.e., absolute_plot object)
absolute_plot <- ggplot(data= absolute_df, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  scale_fill_manual(values = phylum_colors_abs) +
  labs(title = "Taxonomic diversity of absolute abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
absolute_plot
ggsave(absolute_plot, filename = here("absolute_plot.jpg"), width = 10, height = 5)

# Create the figure for the relative abundances
percentages_df$Phylum <- as.factor(percentages_df$Phylum)
phylum_colors_rel<- colorRampPalette(brewer.pal(8,"Dark2")) (length(levels(percentages_df$Phylum)))
relative_plot <- ggplot(data=percentages_df, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  scale_fill_manual(values = phylum_colors_rel) +
  labs(title = "Taxonomic diversity of relative abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
relative_plot
ggsave(relative_plot, filename = here("relative_plot.jpg"), width = 10, height = 5)


# Transformation of the data; Manipulation of the information; and plotting ----
# Subset to the most dominant phylum Bacteroidetes 
bacteriod <- subset_taxa(fmt_metagenomes, Phylum == "Bacteroidetes")
# Transform to relative abundance
bacteriod_percentages <- transform_sample_counts(bacteriod, function(x) x * 100 / sum(x))
# Agglomerate at genus level
bacteriod_glom <- tax_glom(bacteriod_percentages, taxrank = "Genus")
# Convert to data frame for ggplot
bacteriod_df <- psmelt(bacteriod_glom)
# Rename low-abundance genera
bacteriod_df$Genus[bacteriod_df$Abundance < 10] <- "Genera < 10.0 abund."
bacteriod_df$Genus <- as.factor(bacteriod_df$Genus)
# Generate color palette
genus_colors_bacteriod <- colorRampPalette(brewer.pal(8, "Dark2"))(length(levels(bacteriod_df$Genus)))
# Plot
plot_bacteriod <- ggplot(data = bacteriod_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = genus_colors_bacteriod) +
  labs(title = "Diversity of Bacteroidetes at genus level inside the samples") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
plot_bacteriod
ggsave(plot_bacteriod, filename = here("plot_bacteriod.jpg"), width = 10, height = 5)

# Functional analysis visualization ----
# List all .tabular annotation files in the data/ folder
files <- list.files("data", pattern = "annotations\\.tabular$", full.names = TRUE)

# Step 2: Extract accession numbers from filenames
accessions <- str_extract(files, "SRS[0-9]+")

# Step 3: Create your group mapping manually
group_mapping <- tibble(
  accession = c(
    "SRS24541665", "SRS24541666", "SRS24541667", "SRS24541668",  # Responder
    "SRS24541674", "SRS24541675", "SRS24541676", "SRS24541677",  # Non_responder
    "SRS24541679", "SRS24541680", "SRS24541681", "SRS24541682"   # Donor
  ),
  group = c(
    rep("Responder", 4),
    rep("Non_responder", 4),
    rep("Donor", 4)
  ))

# Merge file list with group mapping. This gives you a metadata table with three columns: file, accession, and group.
metadata <- tibble(
  file = files,
  accession = accessions
) %>%
  left_join(group_mapping, by = "accession")

# Read in and tag each annotation file
all_annotations <- lapply(1:nrow(metadata), function(i) {
  df <- read_tsv(metadata$file[i], col_types = cols(.default = "c"))
  df$SampleID <- metadata$accession[i]
  df$Group <- metadata$group[i]
  return(df)
})

# Combine into one data frame
combined_annotations <- bind_rows(all_annotations)

# Check the combined data
glimpse(combined_annotations)
table(combined_annotations$Group)
colnames(combined_annotations) # observe column names

# Visualize functional annotation ----
# 1. Top 10 functional annotations by Group
top_annotation_plot <- 
  combined_annotations %>%
  count(Group, Description) %>%
  group_by(Group) %>%
  slice_max(n, n = 10) %>% 
  ggplot(aes(x = reorder(Description, n), y = n, fill = Group)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Top 10 Functional Annotations by Group", x = "Annotation", y = "Count") +
  theme_minimal()
ggsave(top_annotation_plot, filename = here("top 10 functional annotations.jpg"), width = 10, height = 5, dpi = 300)

# 2. Filter and count KEGG pathways
# Filter out empty or NA entries
kegg_data <- combined_annotations %>%
  filter(!is.na(KEGG_Pathway) & KEGG_Pathway != "-") %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%  # In case multiple pathways are in one cell
  count(Group, KEGG_Pathway, sort = TRUE)

# Visualize top 10 pathways per group
top_kegg_plot <- 
  kegg_data %>%
  group_by(Group) %>%
  slice_max(n, n = 10) %>%
  ggplot(aes(x = reorder(KEGG_Pathway, n), y = n, fill = Group)) +
  geom_col(position = "dodge") +
  labs(title = "Top KEGG Pathways by Group", x = "KEGG Pathway", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave(top_kegg_plot, filename = here("top 10 kegg pathways.jpg"), width = 10, height = 5, dpi = 300)

# 3. Filter and count GO terms
go_data <- combined_annotations %>%
  filter(!is.na(GOs) & GOs != "-") %>%
  separate_rows(GOs, sep = ",") %>%
  count(Group, GOs, sort = TRUE)

# Visualize top 10 GO terms per group
go_data_plot <- 
  go_data %>%
  group_by(Group) %>%
  slice_max(n, n = 10) %>%
  ggplot(aes(x = reorder(GOs, n), y = n, fill = Group)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Top GO Terms by Group", x = "GO Term", y = "Count") +
  theme_minimal()
ggsave(go_data_plot, filename = here("top 10 GO terms.jpg"), width = 10, height = 5, dpi = 300)

# 4. COG functional categories 
cog_functional_plot <- 
  combined_annotations %>%
  filter(!is.na(COG_category) & COG_category != "-") %>%
  separate_rows(COG_category, sep = "") %>%  # Each letter represents a category
  count(Group, COG_category, sort = TRUE) %>%
  ggplot(aes(x = COG_category, y = n, fill = Group)) +
  geom_col(position = "dodge") +
  labs(title = "COG Functional Categories", x = "COG Category", y = "Count") +
  theme_minimal()
ggsave(cog_functional_plot, filename = here("COG functional categories.jpg"), width = 10, height = 5, dpi = 300)

# instead of the script in number 4 above describing COG functional categories with just letters, you can map the meaning of the letters in the plot

# 4. COG functional categories (with functional meanings)
# STEP 1: Create a named vector to map COG letters to functional meanings
cog_map <- c(
  A = "RNA processing",
  B = "Chromatin dynamics",
  C = "Energy production",
  D = "Cell cycle",
  E = "Amino acid metabolism",
  F = "Nucleotide metabolism",
  G = "Carbohydrate metabolism",
  H = "Coenzyme metabolism",
  I = "Lipid metabolism",
  J = "Translation",
  K = "Transcription",
  L = "DNA replication & repair",
  M = "Cell wall biogenesis",
  N = "Cell motility",
  O = "Post-translational mod.",
  P = "Ion transport",
  Q = "Secondary metabolism",
  R = "General prediction",
  S = "Function unknown",
  T = "Signal transduction",
  U = "Intracellular trafficking",
  V = "Defense mechanisms",
  W = "Extracellular structures",
  Y = "Nuclear structure",
  Z = "Cytoskeleton"
)

# STEP 2: Count and clean the COG category data
cog_data <- combined_annotations %>%
  filter(!is.na(COG_category) & COG_category != "-") %>%
  separate_rows(COG_category, sep = "") %>%   # Split multiple letters
  count(Group, COG_category) %>%
  mutate(COG_meaning = cog_map[COG_category]) %>%
  filter(!is.na(COG_meaning))  # Remove unknown letters if any

# STEP 3: Plot the bar chart with COG meanings
 cog_functional_plot <- 
   ggplot(cog_data, aes(x = reorder(COG_meaning, n), y = n, fill = Group)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "COG Functional Categories by Group",
       x = "Functional Category",
       y = "Count") +
  theme_minimal()
ggsave(cog_functional_plot, filename = here("COG functional categories.jpg"), width = 8, height = 5, dpi = 300)

