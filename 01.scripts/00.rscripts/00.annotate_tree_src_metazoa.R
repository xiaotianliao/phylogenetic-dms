#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ggtree")

library(stringr)
library(ggtree)
library(ggplot2)
library(ape)
library(taxizedb)
library(taxize)
library(phytools)
library(rentrez)
library(data.table)

tree <- read.tree("/Users/xl7/Documents/0.Projects/02.phylo-dms/00.inputs/src_metazoa_asr_lgg4/src_metazoa_asr_lgg4.treefile")

midpoint_tree <- midpoint.root(tree, resolve.root = TRUE, outgroup = "NULL")

extract_numbers <- function(strings) {
  # Use sub function with regex to extract the numbers before _0|, _1|, or _2|
  numbers <- sub("_([0-9])\\|.*", "", strings)
  return(numbers)
}
# Extract taxids
taxids <- extract_numbers(midpoint_tree$tip.label)

# Set your ENTREZ API key
api_key <- "018c5dc180b9111f90bccb86bc9c9678f709"  # Replace with your actual API key
set_entrez_key(api_key)

# Define the function to convert taxids to names
convert_taxid_to_name <- function(taxids, chunk_size = 100) {
  # Initialize an empty list to store the results
  names_list <- list()
  
  # Split taxids into chunks
  n <- length(taxids)
  chunks <- split(taxids, ceiling(seq_along(taxids) / chunk_size))
  
  # Process each chunk
  for (chunk in chunks) {
    classifications <- taxize::classification(chunk, db = "ncbi")
    
    for (i in seq_along(chunk)) {
      taxid <- chunk[i]
      classification <- classifications[[i]]
      
      if (!is.null(classification)) {
        family_name <- classification[classification$rank == "order", "name"]
        if (length(family_name) > 0) {
          names_list[[as.character(taxid)]] <- family_name
        } else {
          names_list[[as.character(taxid)]] <- NA
        }
      } else {
        names_list[[as.character(taxid)]] <- NA
      }
    }
  }
  
  return(names_list)
}

# Convert taxids to names
names <- convert_taxid_to_name(as.numeric(taxids), chunk_size = 100)

# Print the results
head(names)

# Convert the list to a data.table
names_dt <- data.table(TaxID = names(names), FamilyName = unlist(names))

# Count the occurrences of each family name
family_counts <- names_dt[, .N, by = FamilyName]
setnames(family_counts, "N", "Count")

family_counts


# Extract taxids from tip labels (adjust this function if needed based on how your tip labels are structured)
extract_numbers <- function(strings) {
  # Assuming tip labels contain taxids embedded before "_0|", "_1|", or "_2|"
  numbers <- sub("_([0-9])\\|.*", "", strings)
  return(numbers)
}

# Extract taxids from tree tip labels
taxids <- extract_numbers(midpoint_tree$tip.label)

# Map order names using the `names` list you generated
order_names <- sapply(taxids, function(x) names[[as.character(x)]])

# Convert the order names to a factor for coloring
order_factors <- as.factor(unlist(order_names))

# Create a data frame for tip labels and their corresponding order
tip_data <- data.frame(label = midpoint_tree$tip.label, order = order_factors)

# Plot the tree with tips colored by order
p <- ggtree(midpoint_tree, aes(color = order)) %<+% tip_data + 
  #geom_tiplab(aes(label = label), hjust = -0.3) + 
  scale_color_manual(values = rainbow(length(unique(order_factors)))) + 
  layout_circular() +theme(legend.position = "none")

# Display the plot
print(p)


nchar('MMTLSVVCPKFTSSDFSFLCTGFLQMGSGMLRAAPHFSPTSALVASGNRGWVVPGLRSEAGSPGARKLGAGELERTRELWE-PEG-LRLKAAPPEL-AATLLLKNAKLGKLLRRLLPLPL-E-CL-LTLEDELAAREARLLWLEGFAGEEALGAPGMRSLMPCLTHHLHSRFASLFLRHSSSRYL-L--PTTAPEDKQGVSPFCSAAPLL-GLLKALAEL--KEAAEETALL-IALLTQESVGPGRYAEAAGEQPGRAERIYGICLDRAPGLREAEW-RLEALLKLALVGRAPGMGDSSSGRPLAPAGEVEGLLFPCVIPTFCKNAQQETRPYPKPKGGLTNDEGQPDYDDEDNWIEGGLDKDSPRGEWSGTRTVQRS-VREGECTR-RSLPGRGAALPQKSRSEGSAAKAKRGGRANPGPAPELLLFHSGGK-GAPTLAEGTPKRAGMETLSLLSSPAKEVLAALAGDSEDPVGCLGLETGPSGEAWRGRGGEGDPATAAAAAAVAAHASPEPQPCGGFGGGGGGGGRSYGAADSKVLWEQGHCGRGASPPADQRYGYSGVSFLTQTFWGAPLAEEGPSLFLKVDLVNRRLLPPRDGLTLLATRPPKEILGPASIPCWGQLLNSLSPVALGRLPQQAARPASSFRCVSRGEIGTHAIQPWGAGEPNYGTGLGGGAAPYTSMSLQCANNLPRDPGDYRAGGALELRGSGLPAGRGAGAQ-PLLEPPGHPVATPALTPGHRGAPPPPLRALFKSENGLLGEHDRSNVRVTAPLPPLQEETPPQPFGGLTSTNALTPPCERPVRRGFGGSMLLRR--TSVLLLKEWEKLMPATSWVLVGLSPSTIPR-RGAAGTPRAGPDAVPAEPRIRLVKGGLEGAPGGSLALGGGGGNVVVGILLERPQFCVLLVGAESGPQL-PPYSLPGPAFLVP-PVTPGFLPIPARFSLRAYLAAGCPSRVTPQAPLSPWELCSRWGHERFVALYTMVPKEPIYIVTDYEARMVAGAVMEIVAAMVVASLALVFVIATEDDDLTFSKGWISEFATPDPLQPSPSVEEPPGDAAAEGDAGKFHARQDTNPAASSSPLLRLGLILNNNCSHPVRRKLNVSPATWFTFRWLQRETQSGDWWEARSLYFEDVSRTQLHRKVGGVPSSGNSGYIPIGEEPFRSNYVAPVDRPLQGGVAMLLGPRRGSIQAFAEEDARYGCCKIAQSDLPPSLEPPWYDFGKIFLSDGNSNFALHCHGRKKSADAERGQLLGLKPGSGPVAHHRPGNPHTHTPHTQGAFLICFFSNEEHDEELLGDGRGESETTKGNYNVSIKVFLDRILRKIQKKYLLSRAPQGGPLEGRRLRGRGGATGFPGAGLESPWGEARWPGPCLVLEDKSCPALEDRSCPAICISCLSRMLFSKDRLPTLSQGGTTAPCLLPAFLFFYSGATLYSLSIRVDWDQKTARGDHIKHYKIRKLDNTGGGYYITTRAQPAPGSCRVP-AFDSVQPELPVRHYMSELLKLSLLCSVIWS-AAKAL-A-THLRARTYARLYLLCF-CLFPHLLKPAPPPLLL-PEEAYSPLLGQSRLLELV-GL-EEEFLPV-LSLSPC-PGVYHPNDPQVGLCYLLTAACTTAMKAQRREGSPQTLGLAVKTKDAWEISRSSIILQRRLGTGCFGDVWMGYNDGLCYLLTKPCPNLKPQTLGLGRDAWEIARETLSLDKKLGQGCFGDVGWLGKADGLCFNLTLVCANLAPQTLGLAKDAWEVAREALLLEKKLSQGVFAEVWCVEFNYGTWNGSTKVANVKTLKPGTMSPKPAFLETAEKGREKPEAQIMCGKPSGKLSWLRHDLKLGVQLYAVVSEKELPIYIVTEFMCVWHLSFFFLSFHDTLLDSQLEEKGLGAVPWGGRQPTPGRVWGRALWFPAVPQVKSHRLVC-ALLSGSMLLAAQRRAENEFLKNAEGQVLKLPQERLLTSLVDMNCNYSAAQKPCWAAVLELGHSRSPCFTKGKQTKEGKRFGWAFICPEAEVDSEPRGRKEPTWGLGPDLRTRIL-LEKQLNTIPVLTLLVTQEGDAVVQGIVLLVSWPIPPPTAWAGLGWGLEGMAYSSRPDSERTAMERMN-YIHRDLRGGERVRGPVAANILVGCWEWLVCKGCLARWRFLRAAQKYLIADFGLPPPARSRSSPAPPLGGGGGQTGAPGPVAGKGGPSGCGSLPCELIEPPTGAASGGSRQRTLPGHQEPLTDEETSVSFSVAVIYNPQQGPALIKQSEFVLGALQMTVIGDLAGAEGMAYIERKNYIHRDLRACSLQLPDLLPPTCTKFPGTAGLISELAKD-WIARVIEDNKWTAPFLLPQVAEGMAYLEAALYGKFTIKSDVWSSSPYGILALTETLCCKIADFGLARIIENEYLAQEELISKGLSWTWIRVYSFEKEIAMERKLVVMNAAPLMYQKFFPYPVFKCWRSLTGLVASIEKRKVVLMMSHRTCNILLFSFAVPPGQLHPSYSVVTKTAHSTSGMLRVSSLYSRGNRAAEVLEQKRALLVEGLIRGGLGHEATDHKLGAPSSSCSLGYHMPCPPLGPPEPGCPAVSRTLYE-RGRLP-A-VEGLGPVDDATFCMEQTWVLAQGNFTPPGR-A-SPA-LLTLALPL--A-L-SLDDLLCNLPEERPTFWRGRGAGGLRKQRLEVGQNLDLGPEYLQSFLEDLVARWPVTQPLCVPA-GLRAVLYFTAWVPNSEVSKPPNQRRTRAAISHLEEPQYESSPRFAHQLAWSSSSQRDAHLEQQQPW--LKTL-AGDQTPPAPSRALGTF')

################################################################################################################
################################################################################################################
################################################################################################################
# Define the function to convert taxids to both species and order names
convert_taxid_to_names <- function(taxids, chunk_size = 100) {
  # Initialize empty lists to store species and order names
  species_list <- list()
  order_list <- list()
  
  # Split taxids into chunks
  chunks <- split(taxids, ceiling(seq_along(taxids) / chunk_size))
  
  # Process each chunk
  for (chunk in chunks) {
    classifications <- taxize::classification(chunk, db = "ncbi")
    
    for (i in seq_along(chunk)) {
      taxid <- chunk[i]
      classification <- classifications[[i]]
      
      if (!is.null(classification)) {
        # Extract species name
        species_name <- classification[classification$rank == "species", "name"]
        if (length(species_name) > 0) {
          species_list[[as.character(taxid)]] <- species_name
        } else {
          species_list[[as.character(taxid)]] <- NA
        }
        
        # Extract order name
        order_name <- classification[classification$rank == "order", "name"]
        if (length(order_name) > 0) {
          order_list[[as.character(taxid)]] <- order_name
        } else {
          order_list[[as.character(taxid)]] <- NA
        }
      } else {
        species_list[[as.character(taxid)]] <- NA
        order_list[[as.character(taxid)]] <- NA
      }
    }
  }
  
  return(list(species = species_list, order = order_list))
}

# Use the function to convert taxids to both species and order names
names_list <- convert_taxid_to_names(as.numeric(taxids), chunk_size = 100)

# Extract species names and replace tip labels
new_tip_labels <- sapply(taxids, function(x) names_list$species[[as.character(x)]])
midpoint_tree$tip.label <- ifelse(!is.na(new_tip_labels), new_tip_labels, midpoint_tree$tip.label)

# Extract order names for coloring
order_names <- sapply(taxids, function(x) names_list$order[[as.character(x)]])
order_factors <- as.factor(unlist(order_names))

# Ensure the `label` column is a character vector
tip_data <- data.frame(label = as.character(midpoint_tree$tip.label), order = order_factors, stringsAsFactors = FALSE)

# Define your color palette
color_palette <- c(
  "Artiodactyla" = "#E41A1C", "Carnivora" = "#377EB8", "Chiroptera" = "#4DAF4A",
  "Cingulata" = "#984EA3", "Dasyuromorphia" = "#FF7F00", "Dermoptera" = "#FFFF33",
  "Didelphimorphia" = "#A65628", "Diprotodontia" = "#F781BF", "Eulipotyphla" = "#999999",
  "Lagomorpha" = "#A6CEE3", "Macroscelidea" = "#1F78B4", "Microbiotheria" = "#B2DF8A",
  "Monotremata" = "#33A02C", "Perissodactyla" = "#FB9A99", "Pholidota" = "#E31A1C",
  "Pilosa" = "#FDBF6F", "Primates" = "#FF7F00", "Proboscidea" = "#CAB2D6",
  "Rodentia" = "#6A3D9A", "Scandentia" = "#FFFF99", "Sirenia" = "#B15928",
  "Tubulidentata" = "#FFED6F", "NA" = "#000000"  # Adjust the NA color if needed
)

# Plot the tree with colored branches and angled, external labels
circular_p <- ggtree(midpoint_tree, branch.length='none', aes(color = order), size = 0.6) %<+% tip_data + 
  geom_tiplab(aes(label = label, angle = angle), geom = "text", offset = 0.1, align = TRUE, size = 3) +
  scale_color_manual(values = color_palette) + 
  layout_circular() +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

# Display the plot
print(circular_p)

# Plot the tree with colored branches and angled, external labels
rec_p <- ggtree(midpoint_tree, branch.length='none', aes(color = order), size = 0.6) %<+% tip_data + 
  geom_tiplab(aes(label = label), geom = "text", offset = 0.1, align = TRUE, size = 3) +
  scale_color_manual(values = color_palette) + 
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

# Display the plot
print(rec_p)
