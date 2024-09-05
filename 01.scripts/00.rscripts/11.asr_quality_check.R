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
input_file <- '/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/05.src_topiary_iqtree/02.out_asr/src_topiary_asr/src_topiary_asr.state'
result <- extract_sequences_and_probabilities(input_file)
sequences <- result$sequences
posterior_probs <- result$posterior_probs
print(sequences)
print(posterior_probs)

posterior_probs$Node16

# Define the bins
bins <- c(0,0.5,0.6,0.7,0.8,0.85, 0.90, 0.95, 1.00)  # Custom bins

# Cut the data into bins
posterior_probs_bins <- cut(posterior_probs$Node34, breaks = bins, right = TRUE, include.lowest = TRUE)

# Create a data frame for ggplot
data_plot <- as.data.frame(table(posterior_probs_bins))

# Plot using ggplot2
p1 <- ggplot(data_plot, aes(x = posterior_probs_bins, y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Posterior Probability", y = "# of Residues", title = "ANC-AS") +
  theme_classic()

# Cut the data into bins
posterior_probs_bins <- cut(posterior_probs$Node35, breaks = bins, right = TRUE, include.lowest = TRUE)

# Create a data frame for ggplot
data_plot <- as.data.frame(table(posterior_probs_bins))

p2 <- ggplot(data_plot, aes(x = posterior_probs_bins, y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Posterior Probability", y = "# of Residues", title = "ANC-S1") +
  theme_classic()


posterior_probs_bins <- cut(posterior_probs$Node37, breaks = bins, right = TRUE, include.lowest = TRUE)

# Create a data frame for ggplot
data_plot <- as.data.frame(table(posterior_probs_bins))

p3 <- ggplot(data_plot, aes(x = posterior_probs_bins, y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Posterior Probability", y = "# of Residues", title = "ANC-S3") +
  theme_classic()

posterior_probs_bins <- cut(posterior_probs$Node39, breaks = bins, right = TRUE, include.lowest = TRUE)

# Create a data frame for ggplot
data_plot <- as.data.frame(table(posterior_probs_bins))

p4 <- ggplot(data_plot, aes(x = posterior_probs_bins, y = Freq)) +
 geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Posterior Probability", y = "# of Residues", title = "ANC-S5") +
  theme_classic()

posterior_probs_bins <- cut(posterior_probs$Node40, breaks = bins, right = TRUE, include.lowest = TRUE)

# Create a data frame for ggplot
data_plot <- as.data.frame(table(posterior_probs_bins))

p5 <- ggplot(data_plot, aes(x = posterior_probs_bins, y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Posterior Probability", y = "# of Residues", title = "ANC-S6") +
  theme_classic()

posterior_probs_bins <- cut(posterior_probs$Node41, breaks = bins, right = TRUE, include.lowest = TRUE)

# Create a data frame for ggplot
data_plot <- as.data.frame(table(posterior_probs_bins))

p6 <- ggplot(data_plot, aes(x = posterior_probs_bins, y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Posterior Probability", y = "# of Residues", title = "ANC-S7") +
  theme_classic()

posterior_probs_bins <- cut(posterior_probs$Node42, breaks = bins, right = TRUE, include.lowest = TRUE)

# Create a data frame for ggplot
data_plot <- as.data.frame(table(posterior_probs_bins))

p7 <- ggplot(data_plot, aes(x = posterior_probs_bins, y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Posterior Probability", y = "# of Residues", title = "ANC-S8") +
  theme_classic()

library(cowplot)

# Combine the plots using plot_grid
combined_plot <- plot_grid(
  p1, p2, p3, p4, p5, p6, p7,    # Specify the plots to be combined
  labels = c("A", "B", "C", "D", "E", "F", "G"),  # Optional labels for each plot
  ncol = 2,    # Number of columns
  align = 'v'  # Align plots vertically
)

# Display the combined plot
print(combined_plot)


























