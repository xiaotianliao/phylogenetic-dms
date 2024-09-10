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
###
# Function to extract sequences and posterior probabilities
extract_sequences_and_probabilities <- function(file_path) {
  # Read the file into a data table
  data <- fread(file_path, header = TRUE, sep = "\t")
  
  # Initialize an empty list to store sequences
  sequences <- list()
  
  # Initialize an empty list to store posterior probabilities
  posterior_probs <- list()
  
  # Loop through each row of the data table
  for (i in 1:nrow(data)) {
    node <- data$Node[i]
    site <- data$Site[i]
    state <- data$State[i]
    
    # Identify the posterior probability for the most likely state
    prob_column <- paste0("p_", state)
    
    # Check if the column exists and has data
    if (prob_column %in% names(data) && !is.na(data[[prob_column]][i])) {
      prob <- data[[prob_column]][i]
    } else {
      prob <- NA  # Use NA to indicate missing probability
    }
    
    # Initialize the sequence and probability list for a new node
    if (!node %in% names(sequences)) {
      sequences[[node]] <- character(max(data$Site))  # Preallocate character vector for sequence
      posterior_probs[[node]] <- numeric(max(data$Site))  # Preallocate numeric vector for probabilities
    }
    
    # Assign the state to the correct site
    sequences[[node]][site] <- state
    
    # Assign the posterior probability to the correct site
    posterior_probs[[node]][site] <- prob
  }
  
  # Convert the lists of states to strings
  for (node in names(sequences)) {
    sequences[[node]] <- paste(sequences[[node]], collapse = "")
  }
  
  return(list(sequences = sequences, posterior_probs = posterior_probs))
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
# Usage example
input_file <- '/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/04.dlg4_topiary/02.out_asr/dlg_topiary_asr/dlg_topiary_asr.state'
result <- extract_sequences_and_probabilities(input_file)
sequences <- result$sequences
posterior_probs <- result$posterior_probs
print(sequences)
print(posterior_probs)













