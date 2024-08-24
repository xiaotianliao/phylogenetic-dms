# Load necessary libraries
library(dplyr)

# Read the CSV file
data <- read.csv("/Users/xl7/Documents/0.Projects/02.phylo-dms/02.outputs/asnc_lg_sequence_comparisons.csv")

# Rename columns for easier reference
colnames(data) <- c("transition", "mutations")

# Display the first few rows of the data
head(data)

data <- data %>%
  mutate(
    start_node = as.numeric(str_extract(transition, "(?<=Node)\\d+(?= to Node)")),
    end_node = as.numeric(str_extract(transition, "(?<=to Node)\\d+"))
  )

head(data)
nrow(data)

# Filter transitions where the difference between start and end nodes is exactly 1
filtered_data <- data %>%
  filter(abs(end_node - start_node) == 1)

# Count the number of mutations for each transition
filtered_data <- filtered_data %>%
  mutate(mutation_count = str_count(mutations, "\\|") + 1)


head(filtered_data)
nrow(filtered_data)

# Save the filtered data to a new CSV file
#write.csv(filtered_data %>% select(transition, mutations, mutation_count), "/Users/xl7/Documents/0.Projects/02.phylo-dms/02.outputs/filtered_sequence_comparisons.csv", row.names = FALSE)

filtered_data$mutation_count

library(viridis)
p1 <- ggplot(filtered_data, aes(x = mutation_count)) +
  geom_histogram(binwidth = 1, fill = viridis(1), color = "black", alpha = 0.1) +
  geom_vline(xintercept = 5, color = "brown3", linetype = "dashed", size = 0.5) +
  ggtitle("AsnC (1757 sequences) - All Transitions") +  
  xlab("Mutation Count") +  
  ylab("Frequency") +  
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) 

p1

# Count the number of transitions with mutation counts below 10
mutation_counts_below_5 <- filtered_data %>%
  filter(mutation_count < 5) %>%
  nrow()

mutation_counts_below_5/nrow(filtered_data)

# Create the histogram
p2 <- ggplot(filtered_data, aes(x = mutation_count)) +
  geom_histogram(binwidth = 1, fill = viridis(1), color = "black", alpha = 0.1) +
  geom_vline(xintercept = 5, color = "brown3", linetype = "dashed", size = 0.5) +
  annotate("text", x = 2,5, y = 380, label = "64.5%", color = "brown3", size = 4) +
  ggtitle("AsnC - Low Mut-order Transitions") +  xlim(0,20) +
  xlab("Mutation Count") +  
  ylab("Frequency") +  
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) 

p2

library(cowplot)
# Combine the plots into a single column layout
combined_plot <- plot_grid(p1, p2, ncol = 1)

# Print the combined plot
print(combined_plot)




