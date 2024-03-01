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
path_apg <- commandArgs(trailingOnly = TRUE)[3]
path_out <- commandArgs(trailingOnly = TRUE)[4]
name <- commandArgs(trailingOnly = TRUE)[5]

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

# Calculating number of tips in the tree
number_tips_in_tree <- length(tree$tip.label)

# Here I need to create some code which calculates the sampling frequency per tip in the tree
# This can be done by finding the genus name of the species in the tree and then finding the number of species in that genus in the wcvp dataset
# Then I can divide the number of tips in the tree by the number of species in the genus in the wcvp dataset

# Creating a data frame of all the species in the tree.
tree_species <- data.frame(species = tree$tip.label)

# Adding a row with the genus name of each species by finding matching row in wcvp dataset
tree_species$genus <- wcvp_accepted_species[match(tree_species$species, wcvp_accepted_species$taxon_name), "genus"]
tree_species$genus <- unlist(tree_species$genus)

# Find out how many species we have sampled in the genera in the tree
number_species_sampled_in_genus <- character(0)
for (genus in unique(tree_species$genus)) {
  number_species_sampled_in_genus <- c(number_species_sampled_in_genus, length(tree_species[which(tree_species$genus == genus), "species"]))
}

# Find out how many species we have in the genera in the wcvp dataset
number_species_in_genus <- character(0)
for (genus in unique(wcvp_accepted_species$genus)) {
  number_species_in_genus <- c(number_species_in_genus, length(wcvp_accepted_species[which(wcvp_accepted_species$genus == genus), "genus"][[1]]))
}

# Create a data frame with genera, number of species sampled, and number of species in the wcvp dataset
genera_sampling_freq <- data.frame(genus = unique(tree_species$genus), number_species_sampled = number_species_sampled_in_genus, number_species_in_wcvp = number_species_in_genus)

# Calculate the sampling frequency for each genus
genera_sampling_freq$sampling_freq <- as.numeric(genera_sampling_freq$number_species_sampled) / as.numeric(genera_sampling_freq$number_species_in_wcvp)

# Now I can create a data frame where I have all the tip labels on one side and then the sampling frequency on the other side
# I have to loop through all genera in tree_species and add the sampling frequency from genera_sampling_freq.
# I can do this by matching the genus name in tree_species with the genus name in genera_sampling_freq and then adding the sampling frequency from genera_sampling_freq to tree_species

# Merge dataframes based on the 'genus' column
merged_data <- merge(tree_species, genera_sampling_freq, by = "genus", all.x = FALSE)

# Now I can select the species name and the sampling fraction columns and then I have the sampling fraction for each species in the tree
sampling_fraction <- merged_data[, c("species", "sampling_freq")]

write.table(sampling_fraction, file = paste0(path_out, name,"_sampling_fraction.txt"), sep = "\t", row.names = FALSE)
