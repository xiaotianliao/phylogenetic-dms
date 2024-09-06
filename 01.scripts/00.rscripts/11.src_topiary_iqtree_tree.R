library(stringr)
library(ggtree)
library(ggplot2)
library(ape)
library(taxizedb)
library(taxize)
library(phytools)
library(rentrez)
library(data.table)
library(dplyr)
library(viridis)


#######################################################################################
###### RECONSTRUCT SEQUENCES FROM .STATE FILE #########################################
#######################################################################################
# Read the FASTA file and extract the sequence headers
fasta_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/02.src_topiary/seed_to_ali/08_alignment_cleaned_trimmed.fasta"
headers <- readLines(fasta_file)[grepl("^>", readLines(fasta_file))]

# Extract unique identifiers and full names
ids_full_names <- data.frame(
  id = str_extract(headers, "(?<=^>)[^|]+"),  # Extracts only the ID without '>'
  full_name = str_extract(headers, "(?<=^>)[^|]+\\|[^|]+\\|[^|]+")  # Extracts the full name
)

# Replace spaces with underscores in the full_name column
ids_full_names$full_name <- gsub(" ", "_", ids_full_names$full_name)

# Display the extracted identifiers and full names for debugging
print("Extracted IDs and full names:")
print(ids_full_names)

# Read the tree file
tree_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/05.src_topiary_iqtree/02.out_asr/src_topiary_asr/src_topiary_asr.treefile"
tree_data <- readLines(tree_file)

# Display the first few lines of the tree data for debugging
print("Original tree data:")
print(head(tree_data))

# Replace the truncated names with the full-length names
for (i in 1:nrow(ids_full_names)) {
  # Ensure pattern matches correctly for replacement
  pattern <- paste0("\\b", ids_full_names$id[i], "\\b\\|[^:]+")  # Adding word boundaries to ensure precise match
  replacement <- ids_full_names$full_name[i]
  
  # Display the pattern and replacement for debugging
  print(paste("Replacing pattern:", pattern, "with:", replacement))
  
  # Perform the replacement
  tree_data <- str_replace_all(tree_data, pattern, replacement)
}

# Display the modified tree data for debugging
print("Modified tree data:")
print(head(tree_data))

output_tree_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/05.src_topiary_iqtree/02.out_asr/src_topiary_asr/updated_tree.treefile"
writeLines(tree_data, output_tree_file)

#######################################################################################
###### RECONSTRUCT SEQUENCES FROM .STATE FILE #########################################
#######################################################################################
tree <- read.tree("/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/05.src_topiary_iqtree/02.out_asr/src_topiary_asr/updated_tree.treefile")
tree

midpoint_tree <- midpoint.root(tree, resolve.root = TRUE, outgroup = "NULL")

# Define a function to extract species names and kinase types from tip labels
extract_species_and_type <- function(strings) {
  # Extract the species name and kinase type from the tip labels
  species_and_type <- sub("^[^|]+\\|([^|]+)\\|([^|]+)$", "\\1|\\2", strings)
  return(species_and_type)
}

# Extract species names and kinase types from the tree's tip labels
species_and_type <- extract_species_and_type(midpoint_tree$tip.label)

# Print the extracted species names and kinase types to check the output
print(species_and_type)

# Replace the tip labels in the tree with the extracted species names and kinase types
midpoint_tree$tip.label <- species_and_type

# Print the first few tip labels to confirm the changes
print(head(midpoint_tree$tip.label))
length(midpoint_tree$tip.label)


# Extract kinase types from the modified tip labels
kinase_types <- sub("^[^|]+\\|([^|]+)$", "\\1", species_and_type)

# Create a data frame for tip labels and their corresponding kinase types
tip_data <- data.frame(label = midpoint_tree$tip.label, kinase_type = kinase_types, stringsAsFactors = FALSE)

# Identify unique kinase types for dynamic coloring
unique_kinase_types <- unique(kinase_types)
print(unique_kinase_types)

# Generate a dynamic color palette based on the number of unique kinase types
color_palette <- setNames(viridis(length(unique_kinase_types)), unique_kinase_types)
print(color_palette)

# Plot the tree with colored branches based on kinase types
p <- ggtree(midpoint_tree, aes(color = kinase_type)) %<+% tip_data + 
  scale_color_manual(values = color_palette) + 
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) + layout_circular() +
  geom_tiplab(aes(label = label, angle = angle), geom = "text", offset = 0.1, align = TRUE, size = 2) +
  labs(color = "Kinase Type")  # Add legend title for clarity

# Display the plot
print(p)


library(ggtree)

# Plot the tree with colored branches based on kinase types
p <- ggtree(midpoint_tree, aes(color = kinase_type)) %<+% tip_data + 
  scale_color_manual(values = color_palette) + 
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  layout_circular() +
  geom_tiplab(aes(label = label, angle = angle), geom = "text", offset = 0.1, align = TRUE, size = 2) +
  labs(color = "Kinase Type")  # Add legend title for clarity

# Label nodes 1 to 10
# Make sure that nodes 1 to 10 are internal nodes in your tree
for (i in 1:100) {
  p <- p + geom_text2(aes(subset = (node == i), label = paste("Node", i)), size = 3, color = "black", hjust = -0.3, vjust = 0.5)
}

# Display the plot
print(p)


rec_p <- gNode35rec_p <- ggtree(midpoint_tree, aes(color = kinase_type)) %<+% tip_data + 
  scale_color_manual(values = color_palette) + 
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  geom_tiplab(aes(label = label), geom = "text", offset = 0.1, align = TRUE, size = 2) +
  labs(color = "Kinase Type")  # Add legend title for clarity

rec_p

find_tec_node <- function(tree, tip_labels) {
  # Find nodes that have "TEC" kinase type, could be customized based on the tree
  tec_tips <- grep("TEC$", tip_labels)  # Find all TEC tips
  mrca <- getMRCA(tree, tec_tips)       # Get the most recent common ancestor of all TEC tips
  return(mrca)
}

# Get the node for TEC clade
tec_node <- find_tec_node(midpoint_tree, midpoint_tree$tip.label)

# Check the identified node
print(tec_node)



# Plot the tree and collapse the TEC clade
rec_p <- ggtree(midpoint_tree, aes(color = kinase_type)) %<+% tip_data + 
  scale_color_manual(values = color_palette) + 
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  geom_tiplab(aes(label = label), size = 2, align = TRUE, offset = 0.1) +
  labs(color = "Kinase Type")  # Add legend title for clarity

# Collapse the TEC clade
rec_p <- rec_p + geom_cladelabel(node = tec_node, label = "TEC clade", offset = 0.2, barsize = 2, fontsize = 3, hjust = 0.5)

# Display the plot with collapsed TEC clade
print(rec_p)

##########
##########
##########
tree <- read.tree("/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/05.src_topiary_iqtree/02.out_asr/src_topiary_asr/updated_tree.treefile")
tree

# Define a function to extract species names and kinase types from tip labels
extract_species_and_type <- function(strings) {
  # Extract the species name and kinase type from the tip labels
  species_and_type <- sub("^[^|]+\\|([^|]+)\\|([^|]+)$", "\\1|\\2", strings)
  return(species_and_type)
}

# Extract species names and kinase types from the tree's tip labels
species_and_type <- extract_species_and_type(tree$tip.label)

# Print the extracted species names and kinase types to check the output
print(species_and_type)

# Replace the tip labels in the tree with the extracted species names and kinase types
tree$tip.label <- species_and_type

# Print the first few tip labels to confirm the changes
print(head(tree$tip.label))
length(tree$tip.label)


# Extract kinase types from the modified tip labels
kinase_types <- sub("^[^|]+\\|([^|]+)$", "\\1", species_and_type)

# Create a data frame for tip labels and their corresponding kinase types
tip_data <- data.frame(label = tree$tip.label, kinase_type = kinase_types, stringsAsFactors = FALSE)

# Identify unique kinase types for dynamic coloring
unique_kinase_types <- unique(kinase_types)
print(unique_kinase_types)

# Generate a dynamic color palette based on the number of unique kinase types
color_palette <- setNames(viridis(length(unique_kinase_types)), unique_kinase_types)
print(color_palette)

p_daylight <- ggtree(tree, aes(color = kinase_type),layout="daylight") %<+% tip_data + 
  scale_color_manual(values = color_palette) + 
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) + 
  #geom_tiplab(aes(label = label), size = 2, align = TRUE, offset = 0.1) +
  labs(color = "Kinase Type")  # Add legend title for clarity

p_daylight


##########
##########
##########

api_key <- "018c5dc180b9111f90bccb86bc9c9678f709"  
set_entrez_key(api_key)

# Extract species names from tip labels
extract_species_name <- function(label) {
  sub("\\|.*$", "", label)  # Extract part before the first "|"
}

species_names <- sapply(midpoint_tree$tip.label, extract_species_name)
species_names

# Extract species names from tip labels
extract_species_name <- function(label) {
  sub("\\|.*$", "", label)  # Extract part before the first "|"
}

species_names <- sapply(midpoint_tree$tip.label, extract_species_name)
species_names
# Use taxize to get order names
convert_species_to_order <- function(species_names, chunk_size = 100) {
  order_list <- list()
  
  # Split species names into chunks
  chunks <- split(species_names, ceiling(seq_along(species_names) / chunk_size))
  
  for (chunk in chunks) {
    classifications <- taxize::classification(chunk, db = "ncbi")
    
    for (i in seq_along(chunk)) {
      species <- chunk[i]
      classification <- classifications[[i]]
      
      if (!is.null(classification)) {
        # Extract order name
        order_name <- classification[classification$rank == "class", "name"]
        if (length(order_name) > 0) {
          order_list[[species]] <- order_name
        } else {
          order_list[[species]] <- NA
        }
      } else {
        order_list[[species]] <- NA
      }
    }
  }
  
  return(order_list)
}

# Get order names for the species
order_names <- convert_species_to_order(species_names)
order_names 

# Ensure the number of order names matches the number of tip labels
# Initialize a vector to store order information
order_vector <- vector("character", length(midpoint_tree$tip.label))
names(order_vector) <- species_names  # Assign species names as names for matching

# Fill order_vector with order names, using NA where order is not found
for (i in seq_along(species_names)) {
  species <- species_names[i]
  if (!is.null(order_names[[species]])) {
    order_vector[i] <- order_names[[species]]
  } else {
    order_vector[i] <- NA
  }
}

# Convert to factors, ensuring the same length as tip labels
order_factors <- as.factor(order_vector)
unique(order_factors)
# Create a data frame with tip labels and their corresponding order
tip_data <- data.frame(label = midpoint_tree$tip.label, order = order_factors, stringsAsFactors = FALSE)

# Check the length to ensure it matches
length(midpoint_tree$tip.label) == length(order_factors)  # This should be TRUE now

# Define a custom color palette for the unique orders
# Define a custom darker color palette for the unique orders
custom_colors <- c(
  # Birds
  "Aves" = "#FF7F0E",             # Bright Orange (Birds)
  
  # Worms
  "Polychaeta" = "#9467BD",       # Soft Violet (Bristle Worms)
  "Pilidiophora" = "#9467BD",     # Soft Violet (Ribbon Worms)
  "Enteropneusta" = "#9467BD",    # Soft Violet (Acorn Worms)
  
  # Insects & Arthropods
  "Insecta" = "#8C564B",          # Deep Burgundy (Insects)
  "Arachnida" = "#8C564B",        # Deep Burgundy (Spiders, Scorpions)
  "Pycnogonida" = "#8C564B",      # Deep Burgundy (Sea Spiders)
  "Thecostraca" = "#8C564B",      # Deep Burgundy (Barnacles)
  "Branchiopoda" = "#8C564B",     # Deep Burgundy (Water Fleas)
  
  # Fish
  "Myxini" = "#1F77B4",           # Soft Blue (Hagfish)
  "Actinopteri" = "#1F77B4",      # Soft Blue (Ray-finned Fishes)
  "Chondrichthyes" = "#1F77B4",   # Soft Blue (Cartilaginous Fishes like Sharks)
  "Cladistia" = "#1F77B4",        # Soft Blue (Bichirs and Reedfishes)
  "Hyperoartia" = "#1F77B4",      # Soft Blue (Lampreys)
  "Leptocardii" = "#1F77B4",      # Gray (Lancelets)
  
  # Sponges
  "Demospongiae" = "#17BECF",     # Light Cyan (Common Sponges)
  "Calcarea" = "#17BECF",         # Light Cyan (Calcareous Sponges)
  "Homoscleromorpha" = "#17BECF", # Light Cyan (Simple Sponges)
  
  # Mammals
  "Mammalia" = "#D62728",         # Bright Red (Mammals, including Primates)
  
  # Clams & Mollusks
  "Bivalvia" = "#2CA02C",         # Forest Green (Clams, Oysters)
  "Gastropoda" = "#2CA02C",       # Forest Green (Snails, Slugs)
  "Polyplacophora" = "#2CA02C",   # Forest Green (Chitons)
  
  # Octopuses & Squids
  "Cephalopoda" = "#9467BD",      # Soft Purple-Blue (Octopuses, Squids)
  
  # Corals & Cnidarians
  "Anthozoa" = "#FFBB78",         # Light Orange (Corals, Sea Anemones)
  "Scyphozoa" = "#FFBB78",        # Light Orange (True Jellyfish)
  "Hydrozoa" = "#FFBB78",         # Light Orange (Hydras, Portuguese Man O' War)
  "Lingulata" = "#FFBB78",        # Gray (Brachiopods)
  
  # Other Echinoderms
  "Holothuroidea" = "#2CA02C",    # Forest Green (Sea Cucumbers)
  "Echinoidea" = "#2CA02C",       # Forest Green (Sea Urchins)
  "Asteroidea" = "#2CA02C",       # Forest Green (Sea Stars)
  
  # Reptiles 
  "Lepidosauria" = "#BCBD22",     # Gray (Lizards, Snakes)
  
  # Other Miscellaneous
  "Filasterea" = "#7F7F7F",       # Gray (Single-celled organisms related to animals)

  # Unknown/NA category
  "<NA>" = "#7F7F7F"              # Yellow-Green (Unknown)
)



# Plot the tree with colored tips by species order
p <- ggtree(midpoint_tree, layout = "circular", aes(color = order),branch.length = "none") %<+% tip_data + 
  #scale_color_manual(values = viridis(length(unique(order_factors)))) + 
  scale_color_manual(values = custom_colors) + 
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) + 
  geom_tiplab(aes(label = label), size = 2, align = TRUE, offset = 0.1) +
  labs(color = "Species Class")  # Add legend title for clarity

# Display the plot
p

# Find the nodes for TEC, SRC, and ABL clades
find_clade_node <- function(tree, clade_tips) {
  clade_node <- getMRCA(tree, which(tree$tip.label %in% clade_tips))
  return(clade_node)
}

# Replace with actual tip labels that correspond to TEC, SRC, and ABL
tec_tips <- grep("TEC", midpoint_tree$tip.label, value = TRUE)
src_tips <- grep("SRC", midpoint_tree$tip.label, value = TRUE)
abl_tips <- grep("ABL", midpoint_tree$tip.label, value = TRUE)

grep("elegans", midpoint_tree$tip.label, value = TRUE)

# Get the nodes for each clade
tec_node <- find_clade_node(midpoint_tree, tec_tips)
src_node <- find_clade_node(midpoint_tree, src_tips)
abl_node <- find_clade_node(midpoint_tree, abl_tips)

# Plot the tree with species orders and annotated clades
p <- ggtree(midpoint_tree, layout = "circular", aes(color = order),branch.length = "none") %<+% tip_data + 
  scale_color_manual(values = custom_colors) + 
  theme(
    #legend.text = element_text(size = 10),
    legend.position = "none",
    #legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) + 
  geom_tiplab(aes(label = label), size = 1.5, align = TRUE, offset = 0.1) +
  labs(color = "Species Class") +  # Add legend title for clarity
  
  # Annotate the TEC clade
  geom_cladelabel(node = tec_node, label = "TEC", offset = 0.01, barsize = 1, fontsize = 4, color = "#377eb8") +
  
  # Annotate the SRC clade
  geom_cladelabel(node = src_node, label = "SRC", offset = 0.01, barsize = 1, fontsize = 4, color = "#4daf4a") +
  
  # Annotate the ABL clade
  geom_cladelabel(node = abl_node, label = "ABL", offset = 0.01, barsize = 1, fontsize = 4, color = "#e41a1c")

# Display the plot with clade annotations
print(p)


# Define the categories and corresponding updated colors
categories <- c(
  "Birds", 
  "Worms", 
  "Insects & Arthropods", 
  "Fish", 
  "Sponges", 
  "Mammals", 
  "Clams & Mollusks", 
  "Octopuses & Squids", 
  "Corals & Cnidarians", 
  "Other Echinoderms", 
  "Miscellaneous"
)

colors <- c(
  "#FF7F0E",  # Bright Orange (Birds)
  "#9467BD",  # Soft Violet (Worms)
  "#8C564B",  # Deep Burgundy (Insects & Arthropods)
  "#1F77B4",  # Soft Blue (Fish)
  "#17BECF",  # Light Cyan (Sponges)
  "#D62728",  # Bright Red (Mammals)
  "#2CA02C",  # Forest Green (Clams & Mollusks)
  "#6A3D9A",  # Soft Purple-Blue (Octopuses & Squids)
  "#FFBB78",  # Light Orange (Corals & Cnidarians)
  "#33A02C",  # Green (Other Echinoderms)
  "#7F7F7F"   # Gray (Miscellaneous)
)

# Create an empty plot with no axes
plot(1:10, type = "n", xlab = "", ylab = "", axes = FALSE)

# Add a legend to the plot
legend("center", legend = categories, fill = colors, cex = 0.8, title = "Categories", box.lty = 0)


ggtree(midpoint_tree, layout = "rectangular", aes(color = order),branch.length = "none") %<+% tip_data + 
  scale_color_manual(values = custom_colors) + 
  theme(
    #legend.text = element_text(size = 10),
    legend.position = "none",
    #legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) + 
  geom_tiplab(aes(label = label), size = 1.5, align = TRUE, offset = 0.1) +
  labs(color = "Species Class") +  # Add legend title for clarity
  
  # Annotate the TEC clade
  geom_cladelabel(node = tec_node, label = "TEC", offset = 0.01, barsize = 1, fontsize = 4, color = "#377eb8") +
  
  # Annotate the SRC clade
  geom_cladelabel(node = src_node, label = "SRC", offset = 0.01, barsize = 1, fontsize = 4, color = "#4daf4a") +
  
  # Annotate the ABL clade
  geom_cladelabel(node = abl_node, label = "ABL", offset = 0.01, barsize = 1, fontsize = 4, color = "#e41a1c")




# Create a data frame with tip labels and their corresponding order
tip_data <- data.frame(label = tree$tip.label, order = order_factors, stringsAsFactors = FALSE)

# Check the length to ensure it matches
length(tree$tip.label) == length(order_factors)  # This should be TRUE now


# Display the plot
p_daylight <- ggtree(tree, aes(color = order),layout="daylight") %<+% tip_data + 
  scale_color_manual(values = viridis(length(unique(order_factors)))) + 
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) + 
  #geom_tiplab(aes(label = label), size = 2, align = TRUE, offset = 0.1) +
  labs(color = "Species Class")  # Add legend title for clarity

p_daylight 
p_daylight + 
  geom_tiplab(aes(label = label), size = 1.5, align = TRUE, offset = 0.6, hjust = 0)






