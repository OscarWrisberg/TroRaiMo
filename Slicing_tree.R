# Load the required packages
library(ape)
library(data.table)
library(dplyr)
library(phytools)

# Command line arguments
input_file_tree <- commandArgs(trailingOnly = TRUE)[1]
output <- commandArgs(trailingOnly = TRUE)[2]
input_file_wcvp <- commandArgs(trailingOnly = TRUE)[3]
path_out <- commandArgs(trailingOnly = TRUE)[4]
apg <- commandArgs(trailingOnly = TRUE)[5]

# Read the WCVP names file into a data frame
cat("Opening ", input_file_wcvp, "\n")
wcvp <- readRDS(input_file_wcvp)

# Read the GBMB tree
cat("Opening ", input_file_tree, "\n")
tree <- read.tree(input_file_tree)

tree$tip.label <- gsub("_", " ", tree$tip.label)
tree$tip.label <- gsub('"', '', tree$tip.label)  # nolint

# Get tip names from the tree
tip_names <- tree$tip.label

# Find matching and non-matching tips
matching_tips <- tip_names[tip_names %in% wcvp$taxon_name] # 76935 tips are matching
not_matching_tips <- tip_names[!(tip_names %in% wcvp$taxon_name)] # Only 2939 tips are not matching

# I also need to do something about the tips that are not in the WCVP file but are in the tree.
# In an ideal world I would like to use the taxonomy_matcher to find the correct name for all the tips in the tree
# But that is for a different day.....

# Finding length of matching tips
cat("Length of matching tips ", length(matching_tips), "\n")
cat("Length of non-matching tips ", length(not_matching_tips), "\n")
# Create a data frame with tip names and families

# Writing the non matching tips to a file.
write.table(not_matching_tips, file.path(path_out, "not_matching_tips.txt"), sep = "\t", row.names = FALSE)


#########################################################################################################
find_family <- function(name_list, wcvp) {
  names <- character(0)
  families <- character(0)
  for (i in seq_along(name_list)) {
    #print progress
    if (!i %% 500) cat("Percentage done", format(round((i / length(name_list)) * 100, 2), nsmall = 2), " at ", format(Sys.time(), '%H:%M:%S'), "\n")

    family <- wcvp[wcvp$taxon_name == name_list[i], "family"]
    family <- as.character(family[1])
    names <- c(names, name_list[i])
    families <- c(families, family)
    #cat(name_list[i], " is found in family ", family, " \n")
  }

  cat("Length names ", length(names), "\n ")
  cat("Length families ", length(families), "\n ")

  df_families <- data.frame(name = names, families = families)
  return(df_families)
}

tips_families <- find_family(matching_tips, wcvp)

head(tips_families)
# Write tip families to a text file
write.table(tips_families, "tips_families.txt", sep = "\t", row.names = FALSE) #nolint
##########################################################################################################

# Loading tips families so I dont have to wait so fucking long..
tips_families <- fread("tips_families.txt")

# Find unique families
unique_families <- unique(tips_families$families)
unique_families <- as.character(unique_families)
cat("This is the number of unique_families \n")
cat(nrow(unique_families), "\n ")

# Create a data frame to store the number of tips in each family
df_number_tips <- data.frame(family = character(0), number_tips = numeric(0))
non_mono_family <- character(0)

###########################
# Pruning tree to families#
###########################
cat("Names tips_families \n")
cat(names(tips_families), " \n")
cat(which(names(tips_families) == "name"), "  ")
cat(which(names(tips_families) == "families"), " \n ")

for (i in seq_along(unique_families)) {
  family_subset <- tips_families[tips_families$families == unique_families[i], ]
  tips_family <- as.character(family_subset$name)

  # Selecting the tips NOT in the family
  tips_not_in_family <- tree$tip.label[!which(tree$tip.label %in% tips_family)]

  # Check if the tips which are in the family ACTUALLY form a monophyletic group in the tree
  # If they do not form a monophyletic group then I will not prune the tree
  if (is.monophyletic(tree, tips_family) == FALSE) {
    #cat("The tips in ", unique_families[i], " do not form a monophyletic group in the tree \n")
    non_mono_family <- c(non_mono_family, unique_families[i])
    next
  }

  # Prune the tree so it only contains the tips in the family
  pruned_tree <- drop.tip(tree, tip = tree$tip.label[!tree$tip.label %in% tips_family])

  # Calculate the number of tips in the pruned tree
  number_tips <- length(pruned_tree$tip.label)

  # Append the number of tips and family to the data frame
  df_number_tips <- rbind(df_number_tips, data.frame(family = unique_families[i], number_tips = number_tips))

  # Save the pruned tree to a file
  write.tree(pruned_tree, paste0(path_out, "pruned_tree_family_", unique_families[i], "_GBMB.txt"))
}

# Save the data frame to a text file
cat("Writing out file to ", file.path(path_out, output), "\n")
write.table(df_number_tips, file.path(path_out, "GBMB_sp_per_families.txt"), sep = "\t", row.names = FALSE)


# Now Ill link the families to the orders from the apgweb_parsed.csv file in order to also add the order to the data frame.
# This will then be used to do the exact same pruning as above but with the orders instead of the families.

# Loading the apgweb_parsed.csv file
cat("Loading the apgweb_parsed.csv file \n")
apg <- fread(apg)

# Find order function
cat("Creating the find_order function \n")
find_order <- function(fams, apg) {
  fam_list <- character(0)
  orders <- character(0)

  for (family in unique(fams)) {
    order <- apg[which(apg$Syn_Fam == family), "Clade"] # Finding the order of that family in APG file
    order <- as.character(order[1]) # Selecting the order of the family
    cat("family ", family, "Order", order, "\n") # Printing the family and order
    fam_list <- c(fam_list, family)
    orders <- c(orders, order)
  }

  df_orders <- data.frame(family = fam_list, order = orders)
  return(df_orders)
}

# Running the function
cat("Running the function \n")
family_orders <- find_order(unique_families, apg)

# Merging the tips_families and family_orders data frames
cat("Merging the tips_families and family_orders data frames \n")
tips_family_orders <- merge(tips_families, family_orders, by.x = "families", by.y = "family")


# Find unique orders
cat("Finding unique orders \n")
unique_orders <- unique(tips_family_orders$order)

# Create a dataframe to store the number of tips in each order
cat("Creating a dataframe to store the number of tips in each order \n")
df_number_tips_orders <- data.frame(order = character(0), number_tips = numeric(0))
non_mono_order <- character(0)


#########################
# Pruning tree to orders#
#########################
cat("Pruining tree to orders \n")

for (i in seq_along(unique_orders)){

  tips_family <- as.character(family_subset$name)
  order_subset <- tips_family_orders[tips_family_orders$order == unique_orders[i], ]
  tips_order <- as.character(order_subset$name)

  # Check if the tips in each order form a monophyletic clade
  if (is.monophyletic(tree, tips_order) == FALSE) {
    #cat("The tips in ", unique_orders[i], " do not form a monophyletic group in the tree \n")
    non_mono_order <- c(non_mono_order, unique_orders[i])
    next
  } else {

  cat("The tips in ", unique_orders[i], " form a monophyletic group in the tree \n")
  # Prune the tree so it only contains the tips in the order
  pruned_tree <- drop.tip(tree, tip = tree$tip.label[!tree$tip.label %in% tips_order])

  # Calculate the number of tips in the pruned tree
  number_tips <- length(pruned_tree$tip.label)

  # Append the number of tips and order to the data frame
  df_number_tips_orders <- rbind(df_number_tips_orders, data.frame(order = unique_orders[i], number_tips = number_tips))

  # Save the pruned tree to a file
  write.tree(pruned_tree, paste0(path_out, "pruned_tree_order_", unique_orders[i], "_GBMB.txt"))
  }
}

# Save the data frame to a text file
cat("Writing out file to ", file.path(path_out, output), "\n")
write.table(df_number_tips_orders, file.path(path_out, output), sep = "\t", row.names = FALSE)

# Save the non_monophyletic families and orders to a file.
cat("Saving the non_monophyletic families and orders to a file \n")
write.table(non_mono_order, file.path(path_out, "non_mono_order.txt"), sep = "\t", row.names = FALSE)
write.table(non_mono_family, file.path(path_out, "non_mono_family.txt"), sep = "\t", row.names = FALSE)

