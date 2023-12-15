# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("data.table", "ape", "phytools", "geiger", "castor", "MonoPhy")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


# Setting wd
setwd("/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data")




# Loading tree 






# Loading the 4 test trees
arecales_tree <- read.tree("/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/pruned_tree__order_Arecales_GBMB.txt")
zingiberales_tree <- read.tree("/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/pruned_tree__order_Zingiberales_GBMB.txt")
vitales_tree <- read.tree("/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/pruned_tree__order_Vitales_GBMB.txt")
geraniales_tree <- read.tree("/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/twice_pruned_tree_Geraniales_GBMB.txt")


# Reading distribution datasets
arecales_distribution <- fread("Arecales_distribution_data.txt")
zingiberales_distribution <- fread("Zingiberales_distribution_data.txt")
vitales_distribution <- fread("Vitales_distribution_data.txt")
geraniales_distribution <- fread("Geraniales_distribution_data.txt")


# Calculate the sampling frequency in our trees
arecales_sampling_frequency <- length(arecales_tree$tip.label)/length(arecales_distribution$taxon_name)
cat("Sampling frequency is: ", arecales_sampling_frequency, "with ", length(arecales_tree$tip.label), " tips in the tree and ", length(arecales_distribution$taxon_name), " species in the distribution dataset \n")
cat("Proportion of species in Wet tropics: ", sum(length(arecales_distribution$climate_description[which(arecales_distribution$climate_description == 1)]))/length(arecales_distribution$climate_description), "\n")


zingiberales_sampling_frequency <- length(zingiberales_tree$tip.label)/length(zingiberales_distribution$taxon_name)
cat("Sampling frequency is: ", zingiberales_sampling_frequency, "with ", length(zingiberales_tree$tip.label), " tips in the tree and ", length(zingiberales_distribution$taxon_name), " species in the distribution dataset \n")
cat("Proportion of species in Wet tropics: ", sum(length(zingiberales_distribution$climate_description[which(zingiberales_distribution$climate_description == 1)]))/length(zingiberales_distribution$climate_description), "\n")

vitales_sampling_frequency <- length(vitales_tree$tip.label)/length(vitales_distribution$taxon_name)
cat("Sampling frequency is: ", vitales_sampling_frequency, "with ", length(vitales_tree$tip.label), " tips in the tree and ", length(vitales_distribution$taxon_name), " species in the distribution dataset \n")
cat("Proportion of species in Wet tropics: ", sum(length(vitales_distribution$climate_description[which(vitales_distribution$climate_description == 1)]))/length(vitales_distribution$climate_description), "\n")

geraniales_sampling_frequency <- length(geraniales_tree$tip.label)/length(geraniales_distribution$taxon_name)
cat("Sampling frequency is: ", geraniales_sampling_frequency, "with ", length(geraniales_tree$tip.label), " tips in the tree and ", length(geraniales_distribution$taxon_name), " species in the distribution dataset \n")
cat("Proportion of species in Wet tropics: ", sum(length(geraniales_distribution$climate_description[which(geraniales_distribution$climate_description == 1)]))/length(geraniales_distribution$climate_description), "\n")
