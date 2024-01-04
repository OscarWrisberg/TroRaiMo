#Prepare some paths test script locally
# path_to_tree = "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/orders/pruned_tree_order_Zingiberales_GBMB.tre"
# output_name = "test_output.jld2"
# tips_orders_family = "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/tips_family_orders.txt"
# wcvp = "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/wcvp_names_apg_aligned.rds"  # Read the WCVP names file into a data frame"
# order_in_question = "Zingiberales"
# path_apg = "/home/au543206/GenomeDK/Trf_models/TroRaiMo/apgweb_parsed.csv"

# Fetching arguments
# Command line arguments
path_to_tree <- commandArgs(trailingOnly = TRUE)[1]
wcvp <- commandArgs(trailingOnly = TRUE)[2]
order_in_question <- commandArgs(trailingOnly = TRUE)[3]
path_apg <- commandArgs(trailingOnly = TRUE)[4]
path_out <- commandArgs(trailingOnly = TRUE)[5]

# Rcode needed to run to get the sampling frequency
# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("data.table", "ape", "phytools")
# , "geiger", "castor", "MonoPhy"

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


# Command line arguments
input_file_tree <- file.path(path_to_tree)
input_file_wcvp <- file.path(wcvp)
order_in_question <- order_in_question
apg <- file.path(path_apg)

################################################################################################################################################
######################################################-- Loading files --#######################################################################
################################################################################################################################################

# Load the wcvp dataset
wcvp <- readRDS(input_file_wcvp)

# Load the tree
tree <- read.tree(path_to_tree)

# Load apg
apg <- fread(apg)

# Filter the wcvp dataset to include only rows where taxon_status == "Accepted"
wcvp_accepted <- subset(wcvp, taxon_status == "Accepted")
wcvp_accepted_species <- subset(wcvp_accepted, taxon_rank == "Species")

# Remove "_"
tree$tip.label <- gsub("_", " ", tree$tip.label)
tree$tip.label <- gsub('"', '', tree$tip.label)  # nolint

################################################################################################################################################
############################################-- Finding the order for each family --#############################################################
################################################################################################################################################

# Find unique families
unique_families <- unique(wcvp_accepted_species$family)
unique_families <- as.character(unique_families)

# Create a data frame to store the number of tips in each family
df_number_tips <- data.frame(family = character(0), number_tips = numeric(0))
non_mono_family <- character(0)

####################################################################################
####################  Finding the order for each family  ###########################
####################################################################################

# Find order function
find_order <- function(fams, apg) {
  fam_list <- character(0)
  orders <- character(0)

  for (family in unique(fams)) {
    order <- apg[which(apg$Syn_Fam == family), "Clade"] # Finding the order of that family in APG file
    order <- as.character(order[1]) # Selecting the order of the family

    fam_list <- c(fam_list, family)
    orders <- c(orders, order)
  }

  df_orders <- data.frame(family = fam_list, order = orders)
  return(df_orders)
}

# Running the function
family_orders <- find_order(unique_families, apg)
#length(family_orders$order)

# Merging the wcvp and family_orders data frames
wcvp_accepted_species_orders <- merge(wcvp_accepted_species, family_orders, by.x = "family", by.y = "family")

# Creating a subset of wcvp where we only have the accepted species in the order
wcvp_accepted_species_orders_subset <- wcvp_accepted_species_orders[which(wcvp_accepted_species_orders$order == order_in_question),]

# Calculating how many species should be in the order
number_species_in_order <- length(wcvp_accepted_species_orders_subset$species)

# Calculating number of tips in the tree
number_tips_in_tree <- length(tree$tip.label)

# Calculating the sampling frequency
sampling_freq <- number_tips_in_tree/number_species_in_order

# Here I need to create some code which calculates the sampling frequency per tip in the tree
# This can be done by finding the genus name of the species in the tree and then finding the number of species in that genus in the wcvp dataset
# Then I can divide the number of tips in the tree by the number of species in the genus in the wcvp dataset

# Creating a data frame of all the species in the tree.
tree_species <- data.frame(species = tree$tip.label)

# Adding a row with the genus name of each species by finding matching row in wcvp dataset
tree_species$genus <- wcvp_accepted_species[match(tree_species$species, wcvp_accepted_species$taxon_name), "genus"]
tree_species$genus <- unlist(tree_species$genus)

# Creating a list of all genera in the order.
genera_in_order <- unique(wcvp_accepted_species_orders_subset$genus)

# For each genus in the order, find the number of species in that genus in the wcvp dataset
number_species_in_genus <- character(0)
for (genus1 in genera_in_order) {
  print(genus1)
  print(wcvp_accepted_species[which(wcvp_accepted_species$genus == genus1),"genus"])
  print(length(wcvp_accepted_species[which(wcvp_accepted_species$genus == genus1),"genus"][[1]]))
  number_species_in_genus <- c(number_species_in_genus, length(wcvp_accepted_species[which(wcvp_accepted_species$genus == genus1), "genus"][[1]]))
}
#number_species_in_genus

# Find out how many species we have sampled in the genera in the order
number_species_sampled_in_genus <- character(0)
for (genus in genera_in_order) {
  #print(genus)
  #print(tree_species[which(tree_species$genus == genus), "species"])
  #cat("\n\n\n")
  number_species_sampled_in_genus <- c(number_species_sampled_in_genus, length(tree_species[which(tree_species$genus == genus), "species"]))
}

# Join the number_species_in_genus to the genera_in_order data frame
genera_in_order <- data.frame(genera_in_order, number_species_in_genus, number_species_sampled_in_genus)

# Add a column called sampling_freq to the genera_in_order data frame
genera_in_order$sampling_freq <- as.numeric(genera_in_order$number_species_sampled_in_genus)/as.numeric(genera_in_order$number_species_in_genus)
genera_in_order

# Now I can create a data frame where I have all the tip labels on one side and then the sampling frequency on the other side
# I have to loop through all genera in tree_species and add the sampling frequency from genera_in_order.
# I can do this by matching the genus name in tree_species with the genus name in genera_in_order and then adding the sampling frequency from genera_in_order to tree_species

# Merge dataframes based on the 'genera_in_order' column
merged_data <- merge(tree_species, genera_in_order, by.x = "genus", by.y = "genera_in_order", all.x = FALSE)

# Now I can select the species name and the sampling fraction columns and then I have the sampling fraction for each species in the tree
sampling_fraction <- merged_data[,c("species", "sampling_freq")]

write.table(sampling_fraction, file = paste0(path_out,order_in_question,"_sampling_fraction.txt"), sep = "\t", row.names = FALSE)