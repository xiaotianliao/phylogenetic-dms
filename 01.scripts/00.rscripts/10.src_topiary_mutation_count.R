# Load necessary libraries
library(dplyr)
library(data.table)
library(Biostrings)
library(seqinr)
library(seqinr)
library(pbapply)
library(viridis)
library(cowplot)
library(ggplot2)

#######################################################################################
###### RECONSTRUCT SEQUENCES FROM .STATE FILE #########################################
#######################################################################################

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
#input_file <- '/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/05.src_topiary_iqtree/02.out_asr/src_topiary_asr/src_topiary_asr.state'
#sequences <- extract_sequences(input_file)
#head(sequences)

#output_file <- '/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/05.src_topiary_iqtree/02.out_asr/src_topiary_asr/src_topiary_reconstructed_sequences.fasta'
#save_sequences(sequences, output_file)


#######################################################################################
###### RECONSTRUCT SEQUENCES FROM .STATE FILE #########################################
#######################################################################################
### ALIGN WITH MUSCLE IN ALIVIEW

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

# Function to iteratively compare adjacent nodes with a progress bar
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
  all_differences <- vector("list", n - 1)
  idx <- 1
  
  # Compare each pair of adjacent nodes using a progress bar
  pb <- pbapply::pboptions(type = "timer")
  results <- pbapply::pblapply(1:(n - 1), function(i) {
    node1 <- node_names[i]
    node2 <- node_names[i + 1]
    compare_sequences(sequences, node1, node2)
  }, cl = 1) # Set cl > 1 for parallel processing if desired
  
  for (result in results) {
    if (!is.null(result)) {
      all_differences[[idx]] <- result
      idx <- idx + 1
    }
  }
  
  # Combine all results into a single data frame
  all_differences <- do.call(rbind, all_differences)
  return(all_differences)
}

# Example usage
fasta_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/05.src_topiary_iqtree/02.out_asr/src_topiary_asr/src_topiary_reconstructed_sequences_aligned_muscle.fasta"
all_differences <- compare_all_nodes(fasta_file)

# Print the result
head(all_differences)
nrow(all_differences)


# Specify the output file path
output_file <- "/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/05.src_topiary_iqtree/02.out_asr/src_topiary_asr/sequence_comparisons.csv"

# Save the data frame to a CSV file
write.csv(all_differences, file = output_file, row.names = FALSE)

compare_sequences(sequences, "Node34", "Node35")
compare_sequences(sequences, "Node34", "Node155")
compare_sequences(sequences, "Node35", "Node36")
compare_sequences(sequences, "Node36", "Node37")
compare_sequences(sequences, "Node37", "Node38")
compare_sequences(sequences, "Node38", "Node39")
compare_sequences(sequences, "Node39", "Node40")
compare_sequences(sequences, "Node40", "Node41")
compare_sequences(sequences, "Node41", "Node42")
compare_sequences(sequences, "Node42", "Node43")
compare_sequences(sequences, "Node43", "Node44") # IDENTICAL!! ## TO HUMAN 2SRC (two insertions)
sequences$Node34
sequences$Node35
sequences$Node36
sequences$Node37
sequences$Node38
sequences$Node39
sequences$Node40
sequences$Node41
sequences$Node42
sequences$Node43
sequences$Node44

#######################################################################################
###### RECONSTRUCT SEQUENCES FROM .STATE FILE #########################################
#######################################################################################
# Read the CSV file
data <- read.csv("/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/05.src_topiary_iqtree/02.out_asr/src_topiary_asr/sequence_comparisons.csv")

# Rename columns for easier reference
colnames(data) <- c("transition", "mutations")

# Display the first few rows of the data
head(data)
library(stringr)

data <- data %>%
  mutate(
    start_node = as.numeric(str_extract(transition, "(?<=Node)\\d+(?= to Node)")),
    end_node = as.numeric(str_extract(transition, "(?<=to Node)\\d+"))
  )

# Filter transitions where the difference between start and end nodes is exactly 1
filtered_data <- data %>%
  filter(abs(end_node - start_node) == 1)

# Count the number of mutations for each transition
filtered_data <- filtered_data %>%
  mutate(mutation_count = str_count(mutations, "\\|") + 1)


head(filtered_data)
nrow(filtered_data)
filtered_data$mutation_count

# Count the number of transitions with mutation counts below 10
mutation_counts_below_2 <- filtered_data %>%
  filter(mutation_count < 2) %>%
  nrow()
mutation_counts_below_2
mutation_counts_below_2/nrow(filtered_data)


p1 <- ggplot(filtered_data, aes(x = mutation_count)) +
  geom_histogram(binwidth = 1, fill = viridis(1), color = "black", alpha = 0.1) +
  geom_vline(xintercept = 5, color = "brown3", linetype = "dashed", size = 0.5) +
  ggtitle("LysR (1757 seq - pFAM ) - All Transitions") +  
  xlab("Mutation Count") +  
  ylab("Frequency") +  
  theme_classic() + 
  annotate("text", x = 2,5, y = 100, label = "18.6%", color = "brown3", size = 4) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) 

p1

filtered_data


