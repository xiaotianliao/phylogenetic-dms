

library(data.table)

# Function to extract sequences
extract_sequences <- function(file_path) {
  # Read the file into a data table
  data <- fread(file_path, header = TRUE, sep = "\t")
  
  # Initialize an empty list to store sequences
  sequences <- list()
  
  # Loop through each row of the data table
  for (i in 1:nrow(data)) {
    node <- data$Node[i]
    site <- data$Site[i]
    state <- data$State[i]
    
    # Initialize the sequence list for a new node
    if (!node %in% names(sequences)) {
      sequences[[node]] <- character(max(data$Site))  # Preallocate character vector of the appropriate length
    }
    
    # Assign the state to the correct site
    sequences[[node]][site] <- state
  }
  
  # Convert the lists of states to strings
  for (node in names(sequences)) {
    sequences[[node]] <- paste(sequences[[node]], collapse = "")
  }
  
  return(sequences)
}

# Function to save sequences to a file in FASTA format
save_sequences <- function(sequences, output_file) {
  # Open a connection to the output file
  file_conn <- file(output_file, "w")
  
  # Write each sequence to the file in FASTA format
  for (node in names(sequences)) {
    writeLines(paste0(">", node), file_conn)
    writeLines(sequences[[node]], file_conn)
  }
  
  # Close the file connection
  close(file_conn)
}

# Usage
#input_file <- '/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/asnc_asr_lg_.state'
input_file <- '/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/src_mammalia_asr/src_mammalia_asr.state'
sequences <- extract_sequences(input_file)

head(sequences)

output_file <- '/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/src_mammalia_asr/src_mammalia_reconstructed_sequences.fasta'
save_sequences(sequences, output_file)

##############################################
library(Biostrings)
library(seqinr)
# Read the FASTA file
fasta_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/mafft_web_asr.fasta"
sequences <- read.fasta(file = fasta_file, seqtype = "DNA")
sequences

#sequences <- readDNAStringSet(fasta_file, format="fasta")
numeric_names <- as.numeric(sub("Node", "", names(sequences)))
head(numeric_names)

# Reorder the sequences based on sorted numeric names
sequences <- sequences[order(numeric_names)]
head(sequences)

sequences


###################################################################################################################
library(seqinr)
library(pbapply)

# Function to compare sequences of two nodes and store differences
compare_sequences <- function(sequences, node1, node2) {
  # Extract the sequences for the specified nodes
  node1_sequence <- sequences[[node1]]
  node2_sequence <- sequences[[node2]]
  
  # Convert sequences to character vectors
  node1_vec <- unlist(strsplit(paste(node1_sequence, collapse = ""), ""))
  node2_vec <- unlist(strsplit(paste(node2_sequence, collapse = ""), ""))
  
  # Ensure sequences are of the same length
  if (length(node1_vec) != length(node2_vec)) {
    stop("Sequences are of different lengths and cannot be directly compared.")
  }
  
  # Compare sequences and find differences
  differences <- node1_vec != node2_vec
  diff_positions <- which(differences)
  
  if (length(diff_positions) == 0) {
    return(NULL)
  }
  
  # Create a formatted difference column
  formatted_diffs <- paste0(toupper(node1_vec[diff_positions]), diff_positions, toupper(node2_vec[diff_positions]))
  
  # Return the formatted differences
  data.frame(
    Comparison = paste0(node1, " to ", node2),
    Differences = paste(formatted_diffs, collapse = " | "),
    stringsAsFactors = FALSE
  )
}

# Function to iteratively compare all nodes with a progress bar
compare_all_nodes <- function(fasta_file) {
  # Read the sequences from the FASTA file
  sequences <- read.fasta(file = fasta_file, seqtype = "DNA")
  
  # Extract numeric node names and reorder sequences
  numeric_names <- as.numeric(sub("Node", "", names(sequences)))
  sequences <- sequences[order(numeric_names)]
  
  # Get all node names
  node_names <- names(sequences)
  
  # Get the number of nodes
  n <- length(node_names)
  
  # Initialize an empty list to store results
  all_differences <- vector("list", n * (n - 1) / 2)
  idx <- 1
  
  # Compare each pair of nodes using outer with a progress bar
  pb <- pbapply::pboptions(type = "timer")
  for (i in 1:(n - 1)) {
    results <- pbapply::pblapply((i + 1):n, function(j) {
      node1 <- node_names[i]
      node2 <- node_names[j]
      compare_sequences(sequences, node1, node2)
    }, cl = 1) # Set cl > 1 for parallel processing if desired
    
    for (result in results) {
      if (!is.null(result)) {
        all_differences[[idx]] <- result
        idx <- idx + 1
      }
    }
  }
  
  # Combine all results into a single data frame
  all_differences <- do.call(rbind, all_differences)
  return(all_differences)
}

# Example usage
#fasta_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/mafft_web_asr.fasta"
fasta_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/src_mammalia_asr/src_mafft_web.fasta"
all_differences <- compare_all_nodes(fasta_file)

# Print the result
head(all_differences)

# Assuming all_differences is already created and contains the desired data

# Specify the output file path
output_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/02.outputs/src_sequence_comparisons.csv"

# Save the data frame to a CSV file
write.csv(all_differences, file = output_file, row.names = FALSE)

# Print a confirmation message
cat("Data has been saved to", output_file, "\n")


