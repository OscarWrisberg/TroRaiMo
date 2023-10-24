# Load the required packages
library(ape)
library(data.table)
library(dplyr)

# Command line arguments
input_file_tree <- commandArgs(trailingOnly = TRUE)[1]
output_folder <- commandArgs(trailingOnly = TRUE)[2]
input_file_wcvp <- commandArgs(trailingOnly = TRUE)[3]
path_out <- commandArgs(trailingOnly = TRUE)[4]
apg <- commandArgs(trailingOnly = TRUE)[5]

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

#head(wcvp)

# Create a data frame with tip names and families
find_family <- function(name_list, wcvp) {
  names <- character(0)
  families <- character(0)
  
  for (name in name_list) {
	  family <- wcvp[wcvp$taxon_name == name, "family"] # Carefull Here I am using family and not family.apg
    family <- as.character(family[1])
	  #print(cat("Name ", name, "Family", family, "\n "))
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
#write.table(tips_families, "tips_families.txt", sep = "\t", row.names = FALSE)

# Find unique families
unique_families <- unique(tips_families$families)

# Create a data frame to store the number of tips in each family
df_number_tips <- data.frame(family = character(0), number_tips = numeric(0))

###########################
# Pruning tree to families#
###########################
for (family in unique_families) {
  tips_family <- tips_families[tips_families$family == family, "name"]
  
  # Prune the tree so it only contains the tips in the family
  pruned_tree <- drop.tip(tree, tree$tip.label[!which(tree$tip.label %in% tips_family)])
  
  # Calculate the number of tips in the pruned tree
  number_tips <- length(pruned_tree$tip.label)
  
  # Append the number of tips and family to the data frame
  df_number_tips <- rbind(df_number_tips, data.frame(family = family, number_tips = number_tips))
  
  # Save the pruned tree to a file
  write.tree(pruned_tree, paste0(path_out,"pruned_tree_family_", family, "_GBMB.txt"))
}

# Save the data frame to a text file
cat("Writing out file to ", file.path(path_out,output), "\n")
write.table(df_number_tips, file.path(path_out,"GBMB_sp_per_families.txt"), sep = "\t", row.names = FALSE)


# Now Ill link the families to the orders from the apgweb_parsed.csv file in order to also add the order to the data frame.
# This will then be used to do the exact same pruning as above but with the orders instead of the families.

# Loading the apgweb_parsed.csv file
cat("Loading the apgweb_parsed.csv file \n")
apg <- fread(apg)

# Find order function
find_order <- function(tip_fams, apg) {
  names <- character(0)
  tip_fams$order <- character(0)
  for (family in unique(tip_fams$families)) {
    order <- apg[apg$Syn_Fam == family, "Clade"] # Finding the order of that family in APG file
    order <- as.character(order[1]) # Selecting the order of the family
    print(cat("family ", family, "Order", order, "\n ")) # Printing the family and order
    tip_fams[tip_fams$families == family, "order"] <- order # Adding the order to the tip_fams data frame
  } 

  return(tip_fams)
}

# Running the function
tips_orders <- find_order(tips_families, apg)

# Find unique orders
unique_orders <- unique(apg$Clade)

# Create a dataframe to store the number of tips in each order
df_number_tips_orders <- data.frame(order = character(0), number_tips = numeric(0))

#########################
# Pruning tree to orders#
#########################

for (order in unique_orders){
  tips_order <- tips_orders[tips_orders$order == order, "name"] # Selecting the tips in the order
  
  # Prune the tree so it only contains the tips in the order
  pruned_tree <- drop.tip(tree, tree$tip.label[!which(tree$tip.label %in% tips_order)])
  
  # Calculate the number of tips in the pruned tree
  number_tips <- length(pruned_tree$tip.label)
  
  # Append the number of tips and order to the data frame
  df_number_tips_orders <- rbind(df_number_tips_orders, data.frame(order = order, number_tips = number_tips))
  
  # Save the pruned tree to a file
  write.tree(pruned_tree, paste0(path_out,"pruned_tree__order_", order, "_GBMB.txt"))
}

# Save the data frame to a text file
cat("Writing out file to ", file.path(path_out,output), "\n")
write.table(df_number_tips_orders, file.path(path_out,), sep = "\t", row.names = FALSE)

