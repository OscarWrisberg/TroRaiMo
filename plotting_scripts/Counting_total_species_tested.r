#########################################################################################################################
################################################## Loading packages  ####################################################
#########################################################################################################################

# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("data.table", "ape", "phytools", "geiger", "dplyr", "ggplot2", "viridis","hrbrthemes", "cowplot", "MetBrewer") #

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Setting output folder
output_folder <- "/home/au543206/GenomeDK/Trf_models/workflow/05_figures"

#########################################################################################################################
############################################ Loading ClaDs runs on orders  ##############################################
#########################################################################################################################

# Loading the ClaDs run for the runs that ran on the Orders
# Specify the folder path
folder_path <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/orders"

# Get a list of all Rdata files in the folder
file_list <- list.files(folder_path, pattern = "Clads_output_.*\\.Rdata", full.names = TRUE)

# Create dataframe for tips and tip rate speciation
clads_tip_lambda <- data.frame()

  
# Loop through each file in the folder
for (i in 1:length(file_list)) { # This loop takes atleast 2 hours to run ....

	# Keeping track of the progress
	print(paste("Processing file", i, "of", length(file_list)))

	# Extractinc the file name
	file = file_list[i]

	# Extract the order name from the file name
	order_name <- gsub("Clads_output_(.*)\\.Rdata", "\\1", basename(file))
	
	# Load the Rdata file
	load(file)

	tree <- CladsOutput$tree

	# Append the tip names and tip rate speciation to the dataframe
	clads_tip_lambda <- rbind(clads_tip_lambda, data.frame(order = order_name,
														   tip_label = tree$tip.label,
														   lambda = CladsOutput$lambdatip_map,
														   extinction = CladsOutput$eps_map
														   ))
	
	# print dim of the dataframe to keep track of the progress
	cat("The dataset contains ",dim(clads_tip_lambda)[1], " rows and ", dim(clads_tip_lambda)[2], " columns\n")

	# Remove the CladsOutput object from the environment
	rm(CladsOutput)
}

###############################################################################################################################
############################################--- Clads on Subclades ---#########################################################
###############################################################################################################################
# Loading the ClaDs runs that ran on the sub_order trees
# Specify the folder path
folder_path <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/subset_of_orders"

# Get a list of all Rdata files in the folder
file_list <- list.files(folder_path, pattern = "Clads_output_.*\\.Rdata", full.names = TRUE)

length(file_list)

# Loop through each file in the folder
for (i in 1:length(file_list)) { # This loop takes atleast 2 hours to run ....

	# Keeping track of the progress
	print(paste("Processing file", i, "of", length(file_list)))

	# Extractinc the file name
	file = file_list[i]

	# Extract the order name from the file name
	order_name <- gsub("Clads_output_(.*)\\.Rdata", "\\1", basename(file))
	
	# Load the Rdata file
	load(file)

	tree <- CladsOutput$tree

	# Append the tip names and tip rate speciation to the dataframe
	clads_tip_lambda <- rbind(clads_tip_lambda, data.frame(order = order_name,
														   tip_label = tree$tip.label,
														   lambda = CladsOutput$lambdatip_map,
														   extinction = CladsOutput$eps_map
														   ))
	
	# print dim of the dataframe to keep track of the progress
	cat("The dataset contains ",dim(clads_tip_lambda)[1], " rows and ", dim(clads_tip_lambda)[2], " columns\n")

	# Remove the CladsOutput object from the environment
	rm(CladsOutput)
}

# I need to do something to make sure there are not any duplicate tip labels in my data.
# Ill start by finding a list of all the tip_labels which have duplicates.

# Find the duplicated tip labels
duplicated_tips <- clads_tip_lambda[duplicated(clads_tip_lambda$tip_label), ]
head(duplicated_tips)

# Now I need to loop through all the tip_labels in duplicated tips and find the orders of these tips
# I will then print the orders of the tips to see if they are the same or different
for ( i in 1:nrow(duplicated_tips)) {
	tip <- duplicated_tips$tip_label[i]
	orders <- clads_tip_lambda[clads_tip_lambda$tip_label == tip, "order"]
	print(orders)
}

###############################################################################################################################
############################################--- Loading WCVP data ---#########################################################
###############################################################################################################################

datadir <- "/home/au543206/GenomeDK/Trf_models/data"

# Load the WCVP data
data <- as.data.frame(fread(file.path(datadir,"wcvp_names.csv")))

species <- data[which(data$taxon_rank == "Species"),] # Selecting only species and varieties/subspecies
accepted_species <- species[which(species$taxon_status=="Accepted"),] # Selecting only Accepted species

# Now we find all the genera present in the Clads Tip Lambda data
# We do this by splitting the tip_label column based on " " and selecting the first element
genera <- unique(sapply(strsplit(clads_tip_lambda$tip_label, " "), "[", 1))

# Now we find the number of species in the clads_tip_lambda data
Nr_species_clads <- length(clads_tip_lambda$tip_label)
Nr_species_clads

# Now we find the number of accepted species in these genera in the wcvp.
Nr_species_accepted_accounted_for_in_clads <- length(accepted_species$taxon_name[accepted_species$genus %in% genera])
Nr_species_accepted_accounted_for_in_clads

# Now we find the total number of accepted species in the wcvp
Nr_accepted_species <- length(accepted_species$taxon_name)
Nr_accepted_species

# How big a proportion of the accepted species in the wcvp are accounted for in the clads data
proportion <- Nr_species_accepted_accounted_for_in_clads/Nr_accepted_species
proportion # as of August 12:2024 we account for roughly 55 % of all accepted species in our ClaDs data.

###############################################################################################################################
############################################--- Loading the Smith and Brown tre---#############################################
###############################################################################################################################

# Load the Smith and Brown tree
phylogeny <- read.tree(file.path(datadir, "GBMB_all_accepted_all_species.tre"))

# Get the tip labels of the tree
tip_labels <- phylogeny$tip.label

