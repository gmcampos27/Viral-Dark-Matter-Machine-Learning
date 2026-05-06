#' @title Anellovirus Dark Matter Analysis Pipeline
#' @author Gabriel Montenegro de Campos
#' @description Phylogenetic tree construction and Pairwise Identity analysis for novel Anelloviruses.

# 1. Load Libraries ----
library(ggtree)
library(ape)
library(ggplot2)
library(phangorn)
library(tidyr)
library(dplyr)
library(phytools)
library(treeio)
library(Biostrings)
library(reshape2)
library(stringr)
library(patchwork) # For combining plots

# Set Working Directory
setwd("path/to/tree(files)")

# 2. Phylogenetic Tree Analysis (ggtree) ----

# Load the consensus tree (e.g., from IQ-TREE)
tree <- read.tree("ORF1_CDS_aln.trim.fasta.contree")

# Load associated metadata
metadata <- read.table("anello_metadata.tsv", sep="\t", header = TRUE)

# Clean Study labels for the legend
metadata$Study <- gsub("DarkMatter", "Dark Matter contig", metadata$Study)
metadata$Study <- gsub("NCBI", "NCBI Virus/SCANellome", metadata$Study)

# Rooting and pre-processing
rooted_tree <- midpoint(tree)

# Base Tree Plot (Fan Layout)
p <- ggtree(rooted_tree, layout = "fan", size = 0.5, color = "grey50") %<+% metadata
p$data$ufboot <- as.numeric(as.character(p$data$label))

# Node numbers identification (Run these lines once to find node IDs)
# p + geom_text2(aes(label = node), size = 2) 

# Define nodes for Genus highlighting (Example nodes provided in original script)
generos_nodes <- c(170, 162, 207, 217, 247, 251, 239)
max_x <- max(p$data$x, na.rm = TRUE)

# Customizing the tree
p_tree <- p +
  # Clade Highlights (Genus level)
  geom_hilight(node = 170, fill = "#D4CCE3", alpha = 0.3, max_x + 0.1) + #Beta
  geom_hilight(node = 217, fill = "#D7E4EA", alpha = 0.3, max_x + 0.1) + #gamma
  geom_hilight(node = 162, fill = "#E1EFD6", alpha = 0.4, max_x + 0.1) + #alpha
  geom_hilight(node = 207, fill = "#e0c2a4", alpha = 0.4, max_x + 0.1) + #heto
  geom_hilight(node = 247, fill = "#91A6FF", alpha = 0.4, max_x + 0.1) + #mem
  geom_hilight(node = 251, fill = "#FF88DC", alpha = 0.4, max_x + 0.1) + #lamed
  geom_hilight(node = 239, fill = "#FAFF7F", alpha = 0.4, max_x + 0.1) + #samek
  
  # Tips: Differentiate Dark Matter from References
  
  geom_tippoint(aes(fill = Study), size = 3, shape = 21, color = "black", stroke = 0.2) +
  scale_fill_manual(values = c("Dark Matter contig" = "#FF5154", "NCBI Virus/SCANellome" = "#999999"), name = "Sequences") +
  
  geom_nodepoint(aes(subset = !is.na(ufboot) & ufboot >= 70 & node %in% generos_nodes), 
                 shape = 18, 
                 size = 3) +
  
  geom_treescale(x = 0, y = 0, offset = 1, fontsize = 3, linesize = 1, color = "black") +
  theme(legend.position = "right", 
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 10.5)) +
  
  geom_cladelabel(node = 217, label = "gammatorque",
                  color = "grey10", 
                  barsize = .7, 
                  offset = 0.15,
                  fontsize = 3.5, 
                  fontface = "italic",
                  align = TRUE, 
                  hjust = 1) +
  
  geom_cladelabel(node = 170, label = "betatorque",
                  color = "grey10", 
                  barsize = .7,
                  offset = 0.15,
                  fontface = "italic",
                  fontsize = 3.5,
                  align = TRUE,
                  hjust = -.25) +
  
  geom_cladelabel(node = 162, label = "alphatorque",
                  color = "grey10", 
                  barsize = .7,
                  offset = 0.15,
                  fontface = "italic",
                  fontsize = 3.5,
                  align = TRUE,
                  hjust = -.25) +
  
  geom_cladelabel(node = 207, label = "hetotorque",
                  color = "grey10", 
                  barsize = .7,
                  angle = 0,
                  offset = 0.15,
                  fontface = "italic",
                  fontsize = 3.5,
                  align = TRUE,
                  hjust = 1) +
  
  geom_cladelabel(node = 247, label = "memtorque",
                  color = "grey10", 
                  barsize = .7,
                  offset = 0.15,
                  fontface = "italic",
                  fontsize = 3.5,
                  align = TRUE,
                  hjust = 1) +
  
  geom_cladelabel(node = 251, label = "lamedtorque",
                  color = "grey10", 
                  barsize = .7,
                  offset = 0.15,
                  fontface = "italic",
                  fontsize = 3.5,
                  align = TRUE,
                  hjust = 1) +
  
  geom_cladelabel(node = 239, label = "samektorque",
                  color = "grey10", 
                  barsize = .7,
                  offset = 0.15,
                  fontface = "italic",
                  fontsize = 3.5,
                  align = TRUE,
                  hjust = 1) +
  
  theme(legend.position = "right",
        legend.box = "horizontal", 
        legend.margin = margin(t = 20)) +
  
  
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  
  geom_text2(aes(label = label, subset = !is.na(ufboot) & ufboot >= 70 & node %in% generos_nodes),
             size = 3,
             fontface = "bold",
             nudge_x = -0.05,
             nudge_y = 0.3)

# 3. Pairwise Identity Analysis (Taxonomic Assignment) ----

# Load ORF1 alignment
aln_set <- readDNAStringSet("samek_aln.fasta")
dna_bin <- as.DNAbin(readDNAMultipleAlignment("samek_aln.fasta"))

# Calculate Identity Matrix (1 - raw distance)
dist_mat <- dist.dna(dna_bin, model = "raw", as.matrix = TRUE, pairwise.deletion = TRUE)
identity_mat <- (1 - dist_mat) * 100

# Target Contig for analysis
target_contig <- "NODE_3_length_3165_cov_61917.528821_2" # Simplified for searching

# Transform matrix to long format for ggplot
df_ident <- melt(identity_mat)
colnames(df_ident) <- c("Seq1", "Seq2", "Identity")

# Filter and Clean Data
plot_data <- df_ident %>% 
  filter(grepl(target_contig, Seq1) & Seq1 != Seq2) %>%
  mutate(Accession = str_extract(Seq2, "^[^, ]+"))

# Taxonomy Mapping (Extracting Species from headers)
mapping <- tibble(header = names(aln_set)) %>%
  mutate(Accession = str_extract(header, "^[^, ]+"),
         Species = str_replace_all(str_extract(header, "(?<=SPECIES=)[^,]+"), "_", " "),
         Label = paste0(Accession, " (", Species, ")"))

# Final Data Join and Species Demarcation Logic (ICTV threshold = 69%)
plot_final <- plot_data %>% 
  left_join(mapping, by = "Accession") %>%
  mutate(Status = case_when(
    Identity == max(Identity) ~ "Top Hit",
    Identity >= 69 ~ "Same species (>69%)",
    TRUE ~ "Different species (<69%)"
  ))

# Identity Plot
p_ident <- ggplot(plot_final, aes(x = reorder(Label, Identity), y = Identity)) +
  geom_bar(aes(fill = Status), stat = "identity", width = 0.7, color = "black") +
  scale_fill_manual(values = c("Top Hit"="#EF767A", "Same species (>69%)"="#49BEAA", "Different species (<69%)"="#456990")) +
  geom_text(aes(label = sprintf("%.1f%%", Identity)), hjust = -0.1, size = 3, fontface = "bold") +
  
  # ICTV Threshold Line
  geom_hline(yintercept = 69, color = "firebrick", linetype = "dashed", size = 0.8) +
  annotate("text", x = 1, y = 70, label = "ICTV Threshold (69%)", color = "firebrick", fontface="italic", size=3) +
  
  coord_flip() +
  scale_y_continuous(limits = c(0, 110), expand = c(0,0)) +
  labs(title = "Pairwise Identity Analysis", subtitle = "NODE_3 vs SCANellome Reference Set", y = "Identity (%)", x = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")

# 4. Combine and Export ----

# Combine Tree and Barplot using Patchwork
combined <- (p_tree | p_ident) + 
  plot_annotation(tag_levels = 'A') + 
  plot_layout(widths = c(1.2, 1))

# Save High-Resolution Figure
ggsave("Figure_Final_Phylogeny_Taxonomy.png", combined, width = 400, height = 220, units = "mm", dpi = 300)

print("Analysis complete.!")