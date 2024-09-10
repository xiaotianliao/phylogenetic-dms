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
input_file <- '/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/04.dlg4_topiary/02.out_asr/dlg_topiary_asr/dlg_topiary_asr.state'
sequences <- extract_sequences(input_file)
head(sequences)

output_file <- '/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/04.dlg4_topiary/02.out_asr/dlg_topiary_asr/dlg_topiary_reconstructed_sequences.fasta'
save_sequences(sequences, output_file)

