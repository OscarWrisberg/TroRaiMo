# path_to_tree <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/subset_of_orders/family_phylo_Ebenaceae.tre"
# wcvp <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/wcvp_names_apg_aligned.rds"


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


################################################################################################################################################
######################################################-- Loading files --#######################################################################
################################################################################################################################################

# Load the wcvp dataset
wcvp <- readRDS(input_file_wcvp)

# Load the tree
tree <- read.tree(path_to_tree)

# Filter the wcvp dataset to include only rows where taxon_status == "Accepted"
wcvp_accepted <- subset(wcvp, taxon_status == "Accepted")
wcvp_accepted_species <- subset(wcvp_accepted, taxon_rank == "Species")

# Remove "_"
tree$tip.label <- gsub("_", " ", tree$tip.label)
tree$tip.label <- gsub('"', '', tree$tip.label)  # nolint

################################################################################################################################################
############################################-- Finding the order for each family --#############################################################
################################################################################################################################################

# Create a data frame to store the number of tips in each family
df_number_tips <- data.frame(family = character(0), number_tips = numeric(0))
non_mono_family <- character(0)

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

# Find the families present in the phylogeny
wcvp_accepted_species_tree <- wcvp_accepted_species[which(wcvp_accepted_species$genus %in% tree_species$genus), ]

# Find unique families
unique_families <- unique(wcvp_accepted_species_tree$family)
unique_families <- as.character(unique_families)

# Find out how many species we have sampled in the genera in the tree
number_species_sampled_in_genus <- character(0)
for (genus in unique(tree_species$genus)) {
  number_species_sampled_in_genus <- c(number_species_sampled_in_genus, length(tree_species[which(tree_species$genus == genus), "species"]))
}

# Find out how many species we have in the genera in the wcvp dataset
number_species_in_genus <- tapply(wcvp_accepted_species_tree$genus, wcvp_accepted_species_tree$genus, length)

# Create a data frame with genera, number of species sampled, and number of species in the wcvp dataset
genera_sampling_freq <- data.frame(genus = unique(tree_species$genus), number_species_sampled = number_species_sampled_in_genus, number_species_in_wcvp = number_species_in_genus)

# Calculate the sampling frequency for each genus
genera_sampling_freq$sampling_freq <- as.numeric(genera_sampling_freq$number_species_sampled) / as.numeric(genera_sampling_freq$number_species_in_wcvp)

# Now I can create a data frame where I have all the tip labels on one side and then the sampling frequency on the other side
# I have to loop through all genera in tree_species and add the sampling frequency from genera_sampling_freq.
# I can do this by matching the genus name in tree_species with the genus name in genera_sampling_freq and then adding the sampling frequency from genera_sampling_freq to tree_species

# Merge dataframes based on the 'genus' column
merged_data <- merge(tree_species, genera_sampling_freq, by = "genus", all.x = TRUE)

# Now I can select the species name and the sampling fraction columns and then I have the sampling fraction for each species in the tree
sampling_fraction <- merged_data[, c("species", "sampling_freq")]

write.table(sampling_fraction, file = paste0(path_out, name,"_sampling_fraction.txt"), sep = "\t", row.names = FALSE)
