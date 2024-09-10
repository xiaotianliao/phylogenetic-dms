# Load necessary libraries
library(ape)
library(stringr)

# Read the FASTA file
fasta_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/04.dlg4_topiary/02.out_asr/dlg_topiary_asr/dlg_topiary_reconstructed_sequences.fasta"
fasta_lines <- readLines(fasta_file)
cat("Read FASTA file with", length(fasta_lines), "lines\n")


# Define the input and output files
newick_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/04.dlg4_topiary/02.out_asr/dlg_topiary_asr/updated_tree.treefile"
output_csv_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/04.dlg4_topiary/02.out_asr/dlg_topiary_asr/updated_tree.csv"

# Read the Newick tree
tree <- read.tree(newick_file)
tip_labels <- tree$tip.label
node_labels <- tree$node.label
tip_labels
node_labels
# Extract tree edge data
edge_data <- data.frame(tree$edge)

# Extract node labels (if available)
node_labels <- data.frame(Label = c(tree$tip.label, tree$node.label))
node_labels$ID <- 1:nrow(node_labels)

# Merge node information with edge data
edge_data <- merge(edge_data, node_labels, by.x = "X1", by.y = "ID", all.x = TRUE)
edge_data <- merge(edge_data, node_labels, by.x = "X2", by.y = "ID", all.x = TRUE)

# Add branch lengths (if available)
edge_data$BranchLength <- tree$edge.length
edge_data

# Save the output as CSV
write.csv(edge_data, output_csv_file, row.names = FALSE)

cat("Newick file converted to CSV and saved as:", output_csv_file, "\n")

output_fasta_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/04.dlg4_topiary/02.out_asr/dlg_topiary_asr/updated_reconstructed.fasta"
fasta_lines <- readLines(output_fasta_file)
head(fasta_lines)

# Read the CSV file
csv_data <- read.csv(output_csv_file)
csv_data
# Create a map of original fasta header to new header (NodeID + Species name)
header_map <- setNames(paste0(csv_data$Label.x, "|", csv_data$Label.y), csv_data$Label.x)
header_map

# Read and update the FASTA file
updated_fasta <- c()
for (i in seq(1, length(fasta_lines), by = 2)) {
  # Extract the header (without '>')
  header <- gsub(">", "", fasta_lines[i])
  sequence <- fasta_lines[i + 1]
  
  # Check if the header is present in the CSV mapping
  if (header %in% names(header_map)) {
    new_header <- paste0(">", header_map[[header]])  # Get the new header from the map
  } else {
    new_header <- paste0(">", header)  # Keep the original header if no match found
  }
  
  # Append the updated header and sequence to the new FASTA
  updated_fasta <- c(updated_fasta, new_header, sequence)
}

# Write the updated FASTA file
writeLines(updated_fasta, output_fasta_file)

cat("Updated FASTA file saved as:", output_fasta_file, "\n")

#######################################################################
# Load necessary libraries
library(ape)

# Read the Newick tree file
newick_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/04.dlg4_topiary/02.out_asr/dlg_topiary_asr/updated_tree.treefile"
tree <- read.tree(newick_file)

# Extract tip labels (species names)
tip_labels <- tree$tip.label

# Filter out the tips containing "DLG4"
dlg4_tips <- tip_labels[grepl("DLG4", tip_labels)]
dlg4_tips 
length(dlg4_tips)
length(unique(dlg4_tips))

# Prune the tree to keep only DLG4 tips
pruned_tree <- drop.tip(tree, setdiff(tip_labels, dlg4_tips))
pruned_tree 

# Write the pruned tree back to a Newick file
output_newick <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/04.dlg4_topiary/02.out_asr/dlg_topiary_asr/updated_tree_dlg4.treefile"
write.tree(pruned_tree, file = output_newick)

cat("Pruned DLG4 tree saved to", output_newick, "\n")

# Define a function to extract species names and kinase types from tip labels
# extract_species_and_type <- function(strings) {
#   # Extract the species name and kinase type from the tip labels
#   species_and_type <- sub("^[^|]+\\|([^|]+)\\|([^|]+)$", "\\1|\\2", strings)
#   return(species_and_type)
# }

# Extract species names and kinase types from the tree's tip labels
# species_and_type <- extract_species_and_type(pruned_tree$tip.label)

# Print the extracted species names and kinase types to check the output
# print(species_and_type)

# Replace the tip labels in the tree with the extracted species names and kinase types
# pruned_tree$tip.label <- species_and_type

# Print the first few tip labels to confirm the changes
print(head(pruned_tree$tip.label))
length(pruned_tree$tip.label)
length(unique(pruned_tree$tip.label))

api_key <- "018c5dc180b9111f90bccb86bc9c9678f709"  
set_entrez_key(api_key)

midpoint_tree <- midpoint.root(pruned_tree, resolve.root = TRUE, outgroup = "NULL")


# Extract species names from tip labels
extract_species_name <- function(label) {
  #sub("\\|.*$", "", label)  # Extract part before the first "|"
  sub("^[^|]*\\|([^|]*)\\|.*$", "\\1", label)
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

  # Fish
  "Actinopteri" = "#1F77B4",      # Soft Blue (Ray-finned Fishes)
  "Chondrichthyes" = "#1F77B4",   # Soft Blue (Cartilaginous Fishes like Sharks)
  "Cladistia" = "#1F77B4",        # Soft Blue (Bichirs and Reedfishes)
  
  # Mammals
  "Mammalia" = "#D62728",         # Bright Red (Mammals, including Primates)
  
  # Reptiles 
  "Lepidosauria" = "#BCBD22",     # Gray (Lizards, Snakes)
  
  # Unknown/NA category
  "<NA>" = "#7F7F7F"              # Yellow-Green (Unknown)
)



# Plot the tree with colored tips by species order
p <- ggtree(midpoint_tree, layout = "circular", aes(color = order)) %<+% tip_data + 
  #scale_color_manual(values = viridis(length(unique(order_factors)))) + 
  scale_color_manual(values = custom_colors) + 
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "right"
  ) + 
  geom_tiplab(aes(label = label), size = 2, align = TRUE, offset = 0.01) +
  labs(color = "Species Class")  # Add legend title for clarity

# Display the plot
p

head(csv_data)
head(fasta_file)

# Initialize columns for N-terminal and C-terminal extension
csv_data$N_terminal_extension <- FALSE
csv_data$C_terminal_extension <- FALSE

# Function to extract the complete sequence corresponding to a header
get_sequence_for_header <- function(header) {
  header_line <- grep(paste0(">", header), fasta_lines)
  
  if (length(header_line) > 0) {
    # Collect all sequence lines until the next header or end of file
    seq_lines <- character(0)
    for (j in (header_line + 1):length(fasta_lines)) {
      if (startsWith(fasta_lines[j], ">")) break
      seq_lines <- c(seq_lines, fasta_lines[j])
    }
    # Collapse sequence lines into a single string
    sequence <- paste(seq_lines, collapse = "")
    return(sequence)
  }
  return(NA)
}

# Iterate through each row in the CSV and check for motifs in the FASTA sequence
for (i in 1:nrow(csv_data)) {
  label <- csv_data$Label.x[i]
  sequence <- get_sequence_for_header(label)
  
  if (!is.na(sequence)) {
    # Check for N-terminal motif 'GEEDIPR'
    if (str_detect(sequence, "GEEDIPR")) {
      csv_data$N_terminal_extension[i] <- TRUE
    }
    
    # Check for C-terminal motif 'EYSRFEA'
    if (str_detect(sequence, "EYSRFEA")) {
      csv_data$C_terminal_extension[i] <- TRUE
    }
  }
}

# Save the updated CSV file
write.csv(csv_data, "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/04.dlg4_topiary/02.out_asr/dlg_topiary_asr/updated_tree_with_extension_info.csv", row.names = FALSE)


new_csv_data <- read.csv("/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/04.dlg4_topiary/02.out_asr/dlg_topiary_asr/updated_tree_with_extension_info.csv")
new_csv_data

# Remove rows with NA in Label.y column
new_csv_data <- new_csv_data[!is.na(new_csv_data$Label.y), ]


tree_tips <-pruned_tree$tip.label
head(tree_tips)
length(tree_tips)

new_csv_data_filtered <- new_csv_data[new_csv_data$Label.y %in% tree_tips, ]
nrow(new_csv_data_filtered)

rownames(new_csv_data_filtered) <- new_csv_data_filtered$Label.y

identical(rownames(new_csv_data_filtered), tree$tip.label[tree$tip.label %in% new_csv_data_filtered$Label.y])

# Plot the tree with colored tips by species order
p <- ggtree(midpoint_tree, layout = "circular", aes(color = order)) %<+% tip_data + 
  #scale_color_manual(values = viridis(length(unique(order_factors)))) + 
  scale_color_manual(values = custom_colors) + 
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "right"
  ) + 
  geom_tiplab(aes(label = label), size = 2, align = TRUE, offset = 0.01) +
  #geom_tiplab(aes(label = ifelse(label == "yYnNgIwpxV|Homo_sapiens|DLG4", label, "")), size = 2, align = TRUE, offset = 0.01) + 
  
  labs(color = "Species Class")  # Add legend title for clarity

p
# Display the plot
p1 <- gheatmap(p, new_csv_data_filtered[, "N_terminal_extension", drop = FALSE], 
               offset = 0.01, width = 0.05, colnames_angle = 90, colnames_offset_y = 0.5,
               colnames = FALSE) +
  scale_fill_viridis_d(option="G",  name = "N-terminal Extension")
p1

library(ggnewscale)
p2 <- p1 + new_scale_fill()
gheatmap(p2, new_csv_data_filtered[, "C_terminal_extension", drop = FALSE], 
         offset = 1, width = 0.05, colnames_angle = 90, colnames_offset_y = 0.5,
         colnames = FALSE) +
  scale_fill_viridis_d(option="D",  name = "C-terminal Extension")


######################################################

# Initialize columns for N-terminal and C-terminal extension
csv_data$GEEDIP <- FALSE
csv_data$GEDDIP <- FALSE
csv_data$GDDDIP <- FALSE
csv_data$GEDDYS <- FALSE
csv_data$GDEDIP <- FALSE

# Function to extract the complete sequence corresponding to a header
get_sequence_for_header <- function(header) {
  header_line <- grep(paste0(">", header), fasta_lines)
  
  if (length(header_line) > 0) {
    # Collect all sequence lines until the next header or end of file
    seq_lines <- character(0)
    for (j in (header_line + 1):length(fasta_lines)) {
      if (startsWith(fasta_lines[j], ">")) break
      seq_lines <- c(seq_lines, fasta_lines[j])
    }
    # Collapse sequence lines into a single string
    sequence <- paste(seq_lines, collapse = "")
    return(sequence)
  }
  return(NA)
}

# Iterate through each row in the CSV and check for motifs in the FASTA sequence
for (i in 1:nrow(csv_data)) {
  label <- csv_data$Label.x[i]
  sequence <- get_sequence_for_header(label)
  
  if (!is.na(sequence)) {
    # Check for N-terminal motif 'GEEDIPR'
    if (str_detect(sequence, "GEEDIPR")) {
      csv_data$GEEDIP[i] <- TRUE
    }
    
    # 
    if (str_detect(sequence, "GEDDIP")) {
      csv_data$GEDDIP[i] <- TRUE
    }

    if (str_detect(sequence, "GDDDIP")) {
      csv_data$GDDDIP[i] <- TRUE
    }
    
    if (str_detect(sequence, "GEDDYS")) {
      csv_data$GEDDYS[i] <- TRUE
    }
    
    if (str_detect(sequence, "GDEDIP")) {
      csv_data$GDEDIP[i] <- TRUE
    }
    
  }
}

# Save the updated CSV file
write.csv(csv_data, "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/04.dlg4_topiary/02.out_asr/dlg_topiary_asr/updated_tree_with_extension_info_details.csv", row.names = FALSE)


new_csv_data <- read.csv("/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/04.dlg4_topiary/02.out_asr/dlg_topiary_asr/updated_tree_with_extension_info_details.csv")
new_csv_data

# Remove rows with NA in Label.y column
new_csv_data <- new_csv_data[!is.na(new_csv_data$Label.y), ]


tree_tips <-pruned_tree$tip.label
head(tree_tips)
length(tree_tips)

new_csv_data_filtered <- new_csv_data[new_csv_data$Label.y %in% tree_tips, ]
nrow(new_csv_data_filtered)

rownames(new_csv_data_filtered) <- new_csv_data_filtered$Label.y

identical(rownames(new_csv_data_filtered), tree$tip.label[tree$tip.label %in% new_csv_data_filtered$Label.y])

# Plot the tree with colored tips by species order
p <- ggtree(midpoint_tree, layout = "circular", aes(color = order),branch.length = "none") %<+% tip_data + 
  #scale_color_manual(values = viridis(length(unique(order_factors)))) + 
  scale_color_manual(values = custom_colors) + 
  theme(
    #legend.text = element_text(size = 10),
    #legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "none"
  ) + 
  #geom_tiplab(aes(label = label), size = 2, align = TRUE, offset = 0.01) +
  geom_tiplab(aes(label = ifelse(label == "yYnNgIwpxV|Homo_sapiens|DLG4", label, "")), size = 2, align = TRUE, offset = 0.01) + 
  
  labs(color = "Species Class")  # Add legend title for clarity

p
library(ggnewscale)

new_csv_data_filtered
p1 <- gheatmap(p, new_csv_data_filtered[, "GEEDIP", drop = FALSE], 
               offset = 0.01, width = 0.05, colnames_angle = 90, colnames_offset_y = 0.5,
               colnames = FALSE) +
  scale_fill_viridis_d(option="G",  name = "GEEDIP")
p1


p2 <- p1 + new_scale_fill()
p2 <-gheatmap(p2, new_csv_data_filtered[, "GEDDIP", drop = FALSE], 
         offset = 1, width = 0.05, colnames_angle = 90, colnames_offset_y = 0.5,
         colnames = FALSE) +
  scale_fill_viridis_d(option="D",  name = "GEDDIP")
p2
p3 <- p2 + new_scale_fill()
p3 <- gheatmap(p3, new_csv_data_filtered[, "GDDDIP", drop = FALSE], 
         offset = 2, width = 0.05, colnames_angle = 90, colnames_offset_y = 0.5,
         colnames = FALSE) +
  scale_fill_viridis_d(option="A",  name = "GDDDIP")
p3
p4 <- p3 + new_scale_fill()
p4 <- gheatmap(p4, new_csv_data_filtered[, "GEDDYS", drop = FALSE], 
               offset = 3, width = 0.05, colnames_angle = 90, colnames_offset_y = 0.5,
               colnames = FALSE) +
  scale_fill_viridis_d(option="A",  name = "GEDDYS")
p4

######################################################

api_key <- "018c5dc180b9111f90bccb86bc9c9678f709"  
set_entrez_key(api_key)

newick_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/04.dlg4_topiary/02.out_asr/dlg_topiary_asr/updated_tree.treefile"
# Read the Newick tree
tree <- read.tree(newick_file)

midpoint_tree_full <- midpoint.root(tree, resolve.root = TRUE, outgroup = "NULL")


# Extract species names from tip labels
extract_species_name <- function(label) {
  #sub("\\|.*$", "", label)  # Extract part before the first "|"
  sub("^[^|]*\\|([^|]*)\\|.*$", "\\1", label)
}

species_names <- sapply(midpoint_tree_full$tip.label, extract_species_name)
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
      
      if (is.data.frame(classification)) {  # Check if classification is a data frame
        # Extract the order name if the rank is found
        order_name <- classification[classification$rank == "order", "name"]
        if (length(order_name) > 0) {
          order_list[[species]] <- order_name
        } else {
          order_list[[species]] <- NA
        }
      } else {
        order_list[[species]] <- NA  # Handle cases where classification is NULL
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
order_vector <- vector("character", length(midpoint_tree_full$tip.label))
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
tip_data <- data.frame(label = midpoint_tree_full$tip.label, order = order_factors, stringsAsFactors = FALSE)

# Check the length to ensure it matches
length(midpoint_tree_full$tip.label) == length(order_factors)  # This should be TRUE now

p <- ggtree(midpoint_tree_full, layout = "circular", aes(color = order), branch.length = "none") %<+% tip_data + 
  scale_color_manual(values = rainbow(length(unique(order_factors)))) + 
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "none"
  ) + 
  geom_tiplab(aes(label = ifelse(label %in% c("yYnNgIwpxV|Homo_sapiens|DLG4", "CpQlqtBcsk|Homo_sapiens|DLG1", 
                                              "WynvbQLjFN|Homo_sapiens|DLG3", "ehslqVUJdo|Homo_sapiens|DLG2", 
                                              "GMoDClyZnF|Homo_sapiens|DLG2", "nwpPfhCIKF|Homo_sapiens|DLG3"), 
                                 label, "")), size = 2, align = TRUE, offset = 0.01) + 
  labs(color = "Species Class")


# Display the plot
p



