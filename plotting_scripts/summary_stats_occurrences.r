# First I need to download the occurrences file which contains all the occurrences of the species in the Smith and Brown Phylogeny
library(phytools)
library(dplyr)

# Read the Smith and Brown renamed phylogeny
smb_tree <- read.tree("/home/owrisberg/Trf_models/data/GBMB_pruned.tre")

# remove _ from the tip labels
smb_tree$tip.label <- gsub("_", " ", smb_tree$tip.label)

# Read the occurrences file
occurrence_data <- readRDS("/home/owrisberg/Trf_models/workflow/01_distribution_data/06_Renamed/gbif_renamed.rds")
head(occurrence_data)

dim(occurrence_data)
#

# How many of these occurrences are of a species that is in the Smith and Brown phylogeny?

# First I need to extract the species names from the Smith and Brown phylogeny
smb_species <- smb_tree$tip.label

# Then I need to extract the species names from the occurrence data
occurrence_species <- occurrence_data$species

# Now I can check how many of the species in the occurrence data are in the Smith and Brown phylogeny
sum(occurrence_species %in% smb_species)

# Create a subset of the occurrence data which only contains the species that are in the Smith and Brown phylogeny
smb_occurrences_subset <- occurrence_data[occurrence_species %in% smb_species,]

# Cool so How many species in the smb tree only have 1 occurrence?
# First I need to count the number of occurrences for each species
occurrence_count <- table(occurrence_species_subset)

# Then I can check how many species in the smb tree only have 1 occurrence
sum(smb_species %in% names(occurrence_count)[occurrence_count == 1])

# How many species have less than 10 occurrences?
sum(smb_species %in% names(occurrence_count)[occurrence_count < 10])

# How many species have less than 100 occurrences?
sum(smb_species %in% names(occurrence_count)[occurrence_count < 100])

# how many species have more than 10.000 occurrences?
sum(smb_species %in% names(occurrence_count)[occurrence_count > 10000])

# How many have more than 1000 occurrences?
sum(smb_species %in% names(occurrence_count)[occurrence_count > 1000])

# What is the mean number of occurrences for the species in the smb tree?
mean(occurrence_count)

# what is the median number of occurrences for the species in the smb tree?
median(occurrence_count)

# What is the maximum number of occurrences for the species in the smb tree?
max(occurrence_count)

# What is the minimum number of occurrences for the species in the smb tree?
min(occurrence_count)


occurrence_count


# I need to find a way to sum up the number of species defined as tropical rainforest, non-tropical rainforest, and widespread.
# I probably need to load all the files in workflow_dir+"03_distribution_data/"+Clads_clades[i]+"_states_"+percentages[j]+"_ClaDs.txt"


# Find all files with the pattern *_ClaDs.txt
folder <- "/home/owrisberg/Trf_models/workflow/03_distribution_data/"

files <- list.files(folder, pattern = "*_distribution_data_ClaDs.txt", full.names = TRUE)

# I basically want to read all these files and put them in one large dataframe
# I can do this by using the lapply function

# Read all the files
data <- lapply(files, read.table, header = TRUE, sep = "\t")

# Join all the dataframes together
all_data <- do.call(rbind, data)

# Calculate the number of species deemed as Tropical Rainforest, Non-Tropical Rainforest, and Widespread for the following thresholds: 0.1, 0.25, 0.33, 0.4

# For each of these thrsholds I will create a dataframe.
# this dataframe will have the columns "species", "Tropical Rainforest", "Non-Tropical Rainforest", "Widespread"
# I will then be able to quickly get an overview of the total number of species which are deemed as Tropical Rainforest, Non-Tropical Rainforest, and Widespread for each threshold.

# Define thresholds
thresholds <- c(0.1, 0.25, 0.33, 0.4)

# Initialize an empty list to store dataframes
df_list <- list()

# Loop over each threshold
for (threshold in thresholds) {
  # Create a new dataframe for each threshold
  df <- all_data %>%
    dplyr::mutate(
      Trf = ifelse(proportion_in_tropical_rainforest > threshold, "Present", "Absent"),
      non_trf = ifelse(proportion_outside_tropical_rainforest > threshold, "Present", "Absent")
    ) %>%
    dplyr::select(wcvp_taxon_name, Trf, non_trf)
  
  # Add the dataframe to the list
  df_list[[as.character(threshold)]] <- df
}

# Optionally, you can access each dataframe by threshold
df_list[["0.1"]]   # Dataframe for threshold 0.1
df_list[["0.25"]]  # Dataframe for threshold 0.25
df_list[["0.33"]]  # Dataframe for threshold 0.33
df_list[["0.4"]]   # Dataframe for threshold 0.4

# Calculate the number of species Present in Trf and Absent in non_trf
# For each threshold
for (threshold in thresholds) {
  df <- df_list[[as.character(threshold)]]
  
  # Calculate the number of species Present in Trf and Absent in non_trf
  present_in_trf <- sum(df$Trf == "Present" & df$non_trf == "Absent")
  
  # Calculate the number of species Absent in Trf and Present in non_trf
  absent_in_trf <- sum(df$Trf == "Absent" & df$non_trf == "Present")
  
  # Calculate the number of species Present in both Trf and non_trf
  present_in_both <- sum(df$Trf == "Present" & df$non_trf == "Present")
  
  # Calculate the number of species Absent in both Trf and non_trf
  absent_in_both <- sum(df$Trf == "Absent" & df$non_trf == "Absent")
  
  # Print the results
  cat("Threshold:", threshold, "\n")
  cat("Present in Trf and Absent in non_trf:", present_in_trf, "\n")
  cat("Absent in Trf and Present in non_trf:", absent_in_trf, "\n")
  cat("Present in both Trf and non_trf:", present_in_both, "\n")
  cat("Absent in both Trf and non_trf:", absent_in_both, "\n")
  cat("\n")
}

# Also in the all data dataframe, can we calculate how many species have NA in occurrences_trf
# and how many species have NA in occurrences_non_trf

# Calculate the number of species with NA in occurrences_trf
na_in_occurrences_trf <- sum(is.na(all_data$occurrences_trf))
na_in_occurrences_trf

# Calculate the number of species with NA in occurrences_non_trf
na_in_occurrences_non_trf <- sum(is.na(all_data$occurrences_non_trf))
na_in_occurrences_non_trf
