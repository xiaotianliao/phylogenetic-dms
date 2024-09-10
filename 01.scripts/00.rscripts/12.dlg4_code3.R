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
fasta_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/04.dlg4_topiary/seed_to_ali/07_alignment_trimmed.fasta"
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
tree_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/04.dlg4_topiary/02.out_asr/dlg_topiary_asr/dlg_topiary_asr.treefile"
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

output_tree_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/04.dlg4_topiary/02.out_asr/dlg_topiary_asr/updated_tree.treefile"
writeLines(tree_data, output_tree_file)
