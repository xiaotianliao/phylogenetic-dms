library(ape)

# Read the Newick tree file
newick_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/06.dlg4_blastp/mafft_web_processed/mafft_web_auto_tree.treefile"
tree <- read.tree(newick_file)
tree

midpoint_tree <- midpoint.root(tree, resolve.root = TRUE, outgroup = "NULL")
ggtree(midpoint_tree, layout = "circular")

# Extract species names from tip labels
extract_species_name <- function(label) {
  # Extract species name between underscores
  sub(".*_([A-Za-z]+_[A-Za-z]+)_.*", "\\1", label)
}

species_names <- sapply(midpoint_tree$tip.label, extract_species_name)
species_names

# Use taxize to get order names with error handling
convert_species_to_order <- function(species_names, chunk_size = 100) {
  order_list <- list()
  
  # Split species names into chunks
  chunks <- split(species_names, ceiling(seq_along(species_names) / chunk_size))
  
  for (chunk in chunks) {
    classifications <- taxize::classification(chunk, db = "ncbi")
    
    for (i in seq_along(chunk)) {
      species <- chunk[i]
      classification <- classifications[[i]]
      
      # Check if the result is a data frame and contains the expected columns
      if (is.data.frame(classification) && "rank" %in% colnames(classification) && "name" %in% colnames(classification)) {
        # Extract order name
        order_name <- classification[classification$rank == "class", "name"]
        if (length(order_name) > 0) {
          order_list[[species]] <- order_name
        } else {
          order_list[[species]] <- NA
        }
      } else {
        # Handle case where classification isn't a data frame or doesn't contain rank info
        order_list[[species]] <- NA
      }
    }
  }
  
  return(order_list)
}


# Get order names for the species
order_names <- convert_species_to_order(species_names)
head(order_names)

na_species <- names(order_names[is.na(order_names)])
na_species

# Function to correct species names by capitalizing the first letter
correct_species_name <- function(species) {
  corrected_name <- gsub("(^[a-z])", toupper(substring(species, 1, 1)), species)
  return(corrected_name)
}

convert_species_to_order <- function(species_names, chunk_size = 100, delay = 2) {
  order_list <- list()
  
  # Split species names into chunks
  chunks <- split(species_names, ceiling(seq_along(species_names) / chunk_size))
  
  for (chunk in chunks) {
    classifications <- taxize::classification(chunk, db = "ncbi")
    
    for (i in seq_along(chunk)) {
      species <- chunk[i]
      classification <- classifications[[i]]
      
      if (!is.null(classification) && is.data.frame(classification)) {
        # Extract order name
        order_name <- classification[classification$rank == "class", "name"]
        if (length(order_name) > 0) {
          order_list[[species]] <- order_name
        } else {
          order_list[[species]] <- NA
        }
      } else {
        order_list[[species]] <- NA
        cat("Query failed for:", species, "\n")
      }
    }
    Sys.sleep(delay)  # Add delay between chunks to avoid rate limiting
  }
  
  return(order_list)
}


retry_query <- function(species_name, retries = 3, delay = 1) {
  for (i in seq_len(retries)) {
    result <- taxize::classification(species_name, db = "ncbi")
    if (!is.null(result) && is.data.frame(result)) {
      return(result)
    }
    Sys.sleep(delay)  # Wait before retrying
  }
  return(NA)
}

handle_na_classifications <- function(order_names) {
  na_species <- names(order_names[is.na(order_names)])
  for (species in na_species) {
    corrected_species <- correct_species_name(species)
    order_names[[species]] <- retry_query(corrected_species)
  }
  return(order_names)
}

# Run the function on species with NA classifications
corrected_order_names <- handle_na_classifications(order_names)

# Print out the corrected order names
corrected_order_names

na_species <- names(order_names[is.na(order_names)])
na_species

na_species_new <- names(order_names[is.na(corrected_order_names)])
na_species_new


# Ensure the number of order names matches the number of tip labels
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
length(midpoint_tree$tip.label) == length(order_factors)  # This should be TRUE 

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
  "Lepidosauria" = "#BCBD22",     # Yellow-Green (Lizards, Snakes)
  
  # Frogs, toads
  "Amphibia" = "#9467BD",         # (Frogs, toads, salamanders)
  
  # Unknown/NA category
  "<NA>" = "#7F7F7F"              # Gray (Unknown)
)

# Plot the tree with colored tips by species order
p <- ggtree(midpoint_tree, layout = "circular", aes(color = order), branch.length = "none") %<+% tip_data + 
  scale_color_manual(values = custom_colors) + 
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "right"
  ) + 
  geom_tiplab(aes(label = label), size = 0, align = TRUE, offset = 0.01) +
  labs(color = "Species Class")  # Add legend title for clarity

# Display the plot
p

midpoint_tree$tip.label
# Replace with actual tip labels that correspond to TEC, SRC, and ABL
dlg1_tips <- grep("disks_large_homolog_1|DLG1|Disks_large1|Disks_large_1|discs_large_homolog_1|Disks_large_-like_protein_1|disks_large-like_1|Disks_large-like_protein_1|disks_large_1|disks_large-like_protein_1|discs_large-like_1", midpoint_tree$tip.label, value = TRUE)
dlg2_tips <- grep("disks_large_homolog_2|DLG2|Disks_large2|Disks_large_2|discs_large_homolog_2", midpoint_tree$tip.label, value = TRUE)
dlg3_tips <- grep("disks_large_homolog_3|DLG3|Disks_large3|Disks_large_3|discs_large_homolog_3", midpoint_tree$tip.label, value = TRUE)
dlg4_tips <- grep("disks_large_homolog_4|DLG4|Disks_large4|Disks_large_4|discs_large_homolog_4|discs_large_MAGUK_scaffold_protein_4|Disks_large-like_protein_4|disks_large_4|Disks_large_like_protein_4|disks_large4|_disks_large-like_4", midpoint_tree$tip.label, value = TRUE)

# Find the unmapped tips by subtracting the DLG1, DLG2, DLG3, and DLG4 tips from all the tips
all_tips <- midpoint_tree$tip.label
mapped_tips <- c(dlg1_tips, dlg2_tips, dlg3_tips, dlg4_tips)
na_tips <- setdiff(all_tips, mapped_tips)

# Now you have the na_tips that were not assigned to any of the DLG1, DLG2, DLG3, or DLG4 groups
na_tips
length(na_tips)

length(dlg4_tips)
length(dlg3_tips)
length(dlg2_tips)
length(dlg1_tips)

length(dlg4_tips)+length(dlg3_tips)+length(dlg2_tips)+length(dlg1_tips)+length(na_tips)

annotation_data <- data.frame(
  species = midpoint_tree$tip.label,  # Tip labels from your tree
  clade = ifelse(midpoint_tree$tip.label %in% dlg1_tips, "DLG1",
                 ifelse(midpoint_tree$tip.label %in% dlg2_tips, "DLG2",
                        ifelse(midpoint_tree$tip.label %in% dlg3_tips, "DLG3", "DLG4")))
)

# Convert the clade column into a factor for easy visualization
annotation_data$clade <- as.factor(annotation_data$clade)
nrow(annotation_data)

# Check the number of rows in annotation_data and the tree
print(nrow(annotation_data))  # This should match the number of tip labels in the tree
print(length(tree$tip.label))  # Compare with this to ensure consistency

rownames(annotation_data) <- annotation_data$species
rownames(annotation_data) 


p1 <- gheatmap(p, annotation_data[, "clade", drop=FALSE], 
               offset = 0.01, width = 0.05, colnames_angle = 90, colnames_offset_y = 0.5,
               colnames = FALSE) +
  scale_fill_viridis_d(option="D",  name = "Clade")
p1



# Load the FASTA file
fasta_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/06.dlg4_blastp/mafft_web_processed/mafft_web_auto_aligned.fasta"
fasta_lines <- readLines(fasta_file)


# Extract tree edge data
csv_data <- data.frame(tree.tip.label = tree$tip.label)
csv_data

# Initialize columns for N-terminal and C-terminal extension
csv_data$N_terminal_extension <- FALSE
csv_data$C_terminal_extension <- FALSE

# Function to extract the NCBI ID from the tree tip label
extract_ncbi_id <- function(label) {
  # Extract the NCBI ID (e.g., NP_001308004, XP_010613106)
  match <- regmatches(label, regexpr("NP_[0-9]+|XP_[0-9]+|AAH[0-9]+|KAF[0-9]+", label))
  if (length(match) > 0) {
    return(match)
  } else {
    return(NA)
  }
}

# Modify get_sequence_for_header to match the NCBI ID only (partial match)
get_sequence_for_ncbi_id <- function(ncbi_id) {
  # Perform a partial match using the NCBI ID (e.g., NP_001308004)
  header_line <- grep(paste0(">", ncbi_id), fasta_lines)
  
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

# Test the header extraction and sequence fetching
header <- csv_data$tree.tip.label[10]
ncbi_id <- extract_ncbi_id(header)
sequence <- get_sequence_for_ncbi_id(ncbi_id)

# Output results
print(ncbi_id)
print(sequence)

# Iterate through each row in the CSV, extract NCBI ID, check motifs, and update CSV
for (i in 1:nrow(csv_data)) {
  label <- csv_data$tree.tip.label[i]
  ncbi_id <- extract_ncbi_id(label)  # Extract NCBI ID
  
  if (!is.na(ncbi_id)) {
    sequence <- get_sequence_for_ncbi_id(ncbi_id)  # Get the sequence for the NCBI ID
    
    if (!is.na(sequence)) {
      # Check for N-terminal motif 'GEEDIPR'
      if (grepl("GEEDIP", sequence)) {
        csv_data$N_terminal_extension[i] <- TRUE
      }
      
      # Check for C-terminal motif 'EYSRFEA'
      if (grepl("EYSRFEA", sequence)) {
        csv_data$C_terminal_extension[i] <- TRUE
      }
    }
  }
}

unique(csv_data$N_terminal_extension)
unique(csv_data$C_terminal_extension)
csv_data
nrow(csv_data)
# Save the updated CSV file
output_csv <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/06.dlg4_blastp/mafft_web_processed/extension.csv"
write.csv(csv_data, output_csv, row.names = FALSE)

rownames(csv_data) <- csv_data$tree.tip.label

annotation_data <- data.frame(
  species = midpoint_tree$tip.label,  # Tip labels from your tree
  clade = ifelse(midpoint_tree$tip.label %in% dlg1_tips, "DLG1",
                 ifelse(midpoint_tree$tip.label %in% dlg2_tips, "DLG2",
                        ifelse(midpoint_tree$tip.label %in% dlg3_tips, "DLG3", "DLG4"))),
  N_terminal_extension = csv_data$N_terminal_extension,  # Use appropriate data
  C_terminal_extension = csv_data$C_terminal_extension   # Use appropriate data
)

rownames(annotation_data) <- annotation_data$species
annotation_data$clade <- as.factor(annotation_data$clade)
nrow(annotation_data)

p1 <- gheatmap(p, annotation_data[, "clade", drop=FALSE], 
               offset = 0.01, width = 0.05, colnames_angle = 90, colnames_offset_y = 0.5,
               colnames = FALSE) +
  scale_fill_viridis_d(option="D",  name = "Clade")
p1


p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, annotation_data[, "N_terminal_extension", drop = FALSE], 
               offset = 3, width = 0.05, colnames_angle = 90, colnames_offset_y = 0.5,
               colnames = FALSE) +
  scale_fill_viridis_d(option="G",  name = "N-terminal Extension")
p2

p3 <- p2 + new_scale_fill()
p3 <- gheatmap(p3, annotation_data[, "C_terminal_extension", drop = FALSE], 
               offset = 6, width = 0.05, colnames_angle = 90, colnames_offset_y = 0.5,
               colnames = FALSE) +
  scale_fill_viridis_d(option="G",  name = "C-terminal Extension")
p3




