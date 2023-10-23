# Load the required packages
library(ape)
library(data.table)
library(dplyr)

# Command line arguments
input_file_tree <- commandArgs(trailingOnly = TRUE)[1]
output_file <- commandArgs(trailingOnly = TRUE)[2]
input_file_wcvp <- commandArgs(trailingOnly = TRUE)[3]

# Read the WCVP names file into a data frame
cat("Opening ", input_file_wcvp, "\n")
wcvp <- readRDS(input_file_wcvp)

# Read the GBMB tree
cat("Opening ", input_file_tree, "\n")
tree <- read.tree(input_file_tree)

# Get tip names from the tree
tip_names <- tree$tip.label

# Replace underscores with spaces and remove quotes
tip_names <- gsub("_", " ", tip_names)
tip_names <- gsub('"', '', tip_names)

#cat("Checking the start of the tip_names", tip_names[1:10], "\n")

# Find matching and non-matching tips
matching_tips <- tip_names[tip_names %in% wcvp$taxon_name] # 76935 tips are matching
not_matching_tips <- tip_names[!(tip_names %in% wcvp$taxon_name)] # Only 2939 tips are not matching

# Finding length of matching tips
cat("Length of matching tips ", length(matching_tips), "\n")
cat("Length of non-matching tips ", length(not_matching_tips), "\n")

head(wcvp)

# Create a data frame with tip names and families
find_family <- function(name_list, wcvp) {
  names <- character(0)
  families <- character(0)
  
  for (name in name_list) {
	  family <- wcvp[wcvp$taxon_name == name, "family"] # Carefull Here I am using family and not family.apg
	print(cat("Name ", name, "Family", family, "\n "))
    names <- c(names, name)
    families <- c(families, family)
  }

  print(cat("Length names ", length(names), "\n "))
  print(cat("Length families ", length(families), "\n "))

  df_families <- data.frame(name = names, families = families)
  return(df_families)
}

tips_families <- find_family(matching_tips, wcvp)

# Write tip families to a text file
write.table(tips_families, "tips_families.txt", sep = "\t", row.names = FALSE)

# Find unique families
unique_families <- unique(tips_families$families)

# Create a data frame to store the number of tips in each family
df_number_tips <- data.frame(family = character(0), number_tips = numeric(0))

for (family in unique_families) {
  tips_family <- tips_families[tips_families$family == family, "name"]
  
  # Prune the tree for each family
  pruned_tree <- drop.tip(tree, tips_family)
  
  # Calculate the number of tips in the pruned tree
  number_tips <- length(pruned_tree$tip.label)
  
  # Append the number of tips and family to the data frame
  df_number_tips <- rbind(df_number_tips, data.frame(family = family, number_tips = number_tips))
  
  # Save the pruned tree to a file
  write.tree(pruned_tree, paste0("pruned_tree_", family, "_GBMB.txt"))
}

# Save the data frame to a text file
write.table(df_number_tips, paste0("nr_subtrees_family_slicing_", length(unique_families), ".txt"), sep = "\t", row.names = FALSE)
