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


# I think I need to do this using the Esse trees as I know I finished all runs on these.

# Load the Esse data
# Esse on orders
Esse_order_folder <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/orders"

# Esse on families
Esse_family_folder <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/families"

# Esse on subfamilies
Esse_subfamily_folder <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/subset_of_orders"

# The pattern we want to find in each folder is 
pattern <- "_Esse_tree.tre"

# Get a list of all Rdata files in the folder
file_list_1 <- list.files(Esse_order_folder, pattern = pattern, full.names = TRUE)
file_list_2 <- list.files(Esse_family_folder, pattern = pattern, full.names = TRUE)
file_list_3 <- list.files(Esse_subfamily_folder, pattern = pattern, full.names = TRUE)

file_list <- c(file_list_1, file_list_2, file_list_3)
file_list

# List of orders/fams/subfams in workflow
orders_fams_subfams_list <- list(
    "Acanthaceae_Martyniaceae_Pedaliaceae_1",
    "Acanthaceae_Martyniaceae_Pedaliaceae_2",
    "Acanthaceae_Martyniaceae_Pedaliaceae_3",
    "Acanthaceae_Martyniaceae_Pedaliaceae_4",
    "Acanthaceae_Martyniaceae_Pedaliaceae_5",
    "Amaryllidaceae_1",
    "Amaryllidaceae_2",
    "Anacardiaceae_Burseraceae_Kirkiaceae_1",
    "Anacardiaceae_Burseraceae_Kirkiaceae_2",
    "Apiaceae_1",
    "Apiaceae_2",
    "Apiaceae_3",
    "Apiaceae_4",
    "Apiaceae_5",
    "Apiaceae_6",
    "Apocynaceae_1",
    "Apocynaceae_2",
    "Apocynaceae_3",
    "Apocynaceae_4",
    "Apocynaceae_5",
    "Arecaceae_1",
    "Arecaceae_2",
    "Arecaceae_3",
    "Arecaceae_4",
    "Asparagaceae_1",
    "Asparagaceae_2",
    "Asparagaceae_3",
    "Asteraceae_10",
    "Asteraceae_11",
    "Asteraceae_12",
    "Asteraceae_13",
    "Asteraceae_14",
    "Asteraceae_15",
    "Asteraceae_16",
    "Asteraceae_17",
    "Asteraceae_18",
    "Asteraceae_19",
    "Asteraceae_1",
    "Asteraceae_20",
    "Asteraceae_21",
    "Asteraceae_22",
    "Asteraceae_23",
    "Asteraceae_24",
    "Asteraceae_25",
    "Asteraceae_26",
    "Asteraceae_27",
    "Asteraceae_28",
    "Asteraceae_29",
    "Asteraceae_2",
    "Asteraceae_30",
    "Asteraceae_31",
    "Asteraceae_32",
    "Asteraceae_33",
    "Asteraceae_34",
    "Asteraceae_35",
    "Asteraceae_36",
    "Asteraceae_37",
    "Asteraceae_38",
    "Asteraceae_39",
    "Asteraceae_3",
    "Asteraceae_40",
    "Asteraceae_4",
    "Asteraceae_5",
    "Asteraceae_6",
    "Asteraceae_7",
    "Asteraceae_8",
    "Asteraceae_9",
    "Brassicaceae_1",
    "Brassicaceae_2",
    "Brassicaceae_3",
    "Brassicaceae_4",
    "Bromeliaceae_1",
    "Bromeliaceae_2",
    "Bromeliaceae_3",
    "Bromeliaceae_4",
    "Bromeliaceae_5",
    "Caryophyllaceae_Achatocarpaceae_Amaranthaceae_1",
    "Caryophyllaceae_Achatocarpaceae_Amaranthaceae_2",
    "Caryophyllaceae_Achatocarpaceae_Amaranthaceae_3",
    "Caryophyllaceae_Achatocarpaceae_Amaranthaceae_4",
    "Caryophyllaceae_Achatocarpaceae_Amaranthaceae_5",
    "Caryophyllaceae_Achatocarpaceae_Amaranthaceae_6",
    "Caryophyllaceae_Achatocarpaceae_Amaranthaceae_7",
    "Caryophyllaceae_Achatocarpaceae_Amaranthaceae_8",
    "Cyperaceae_1",
    "Cyperaceae_2",
    "Cyperaceae_3",
    "Cyperaceae_4",
    "Cyperaceae_5",
    "Cyperaceae_6",
    "Cyperaceae_7",
    "Cyperaceae_8",
    "Ericaceae_Clethraceae_Cyrillaceae_1",
    "Ericaceae_Clethraceae_Cyrillaceae_2",
    "Ericaceae_Clethraceae_Cyrillaceae_3",
    "Euphorbiaceae_1",
    "Euphorbiaceae_2",
    "Euphorbiaceae_3",
    "Fabaceae_10",
    "Fabaceae_11",
    "Fabaceae_12",
    "Fabaceae_13",
    "Fabaceae_14",
    "Fabaceae_15",
    "Fabaceae_16",
    "Fabaceae_17",
    "Fabaceae_18",
    "Fabaceae_19",
    "Fabaceae_1",
    "Fabaceae_20",
    "Fabaceae_2",
    "Fabaceae_3",
    "Fabaceae_4",
    "Fabaceae_5",
    "Fabaceae_6",
    "Fabaceae_7",
    "Fabaceae_8",
    "Fabaceae_9",
    "Gesneriaceae_Calceolariaceae_1",
    "Gesneriaceae_Calceolariaceae_2",
    "Gesneriaceae_Calceolariaceae_3",
    "Gesneriaceae_Calceolariaceae_4",
    "Gesneriaceae_Calceolariaceae_5",
    "Gesneriaceae_Calceolariaceae_6",
    "Iridaceae_1",
    "Iridaceae_2",
    "Iridaceae_3",
    "Iridaceae_4",
    "Lamiaceae_10",
    "Lamiaceae_1",
    "Lamiaceae_2",
    "Lamiaceae_3",
    "Lamiaceae_4",
    "Lamiaceae_5",
    "Lamiaceae_6",
    "Lamiaceae_7",
    "Lamiaceae_8",
    "Lamiaceae_9",
    "Lauraceae_1",
    "Lauraceae_2",
    "Lauraceae_3",
    "Melastomataceae_1",
    "Melastomataceae_2",
    "Melastomataceae_3",
    "Melastomataceae_4",
    "Melastomataceae_5",
    "Melastomataceae_6",
    "Melastomataceae_7",
    "Melastomataceae_8",
    "Myrtaceae_1",
    "Myrtaceae_2",
    "Myrtaceae_3",
    "Myrtaceae_4",
    "Myrtaceae_5",
    "Ochnaceae_Clusiaceae_Erythroxylaceae_Podostemaceae_Bonnetiaceae_Rhizophoraceae_Calophyllaceae_Hypericaceae_Ctenolophonaceae_Irvingiaceae_Pandaceae_1",
    "Ochnaceae_Clusiaceae_Erythroxylaceae_Podostemaceae_Bonnetiaceae_Rhizophoraceae_Calophyllaceae_Hypericaceae_Ctenolophonaceae_Irvingiaceae_Pandaceae_2",
    "Ochnaceae_Clusiaceae_Erythroxylaceae_Podostemaceae_Bonnetiaceae_Rhizophoraceae_Calophyllaceae_Hypericaceae_Ctenolophonaceae_Irvingiaceae_Pandaceae_3",
    "Ochnaceae_Clusiaceae_Erythroxylaceae_Podostemaceae_Bonnetiaceae_Rhizophoraceae_Calophyllaceae_Hypericaceae_Ctenolophonaceae_Irvingiaceae_Pandaceae_4",
    "Ochnaceae_Clusiaceae_Erythroxylaceae_Podostemaceae_Bonnetiaceae_Rhizophoraceae_Calophyllaceae_Hypericaceae_Ctenolophonaceae_Irvingiaceae_Pandaceae_5",
    "Orchidaceae_10",
    "Orchidaceae_11",
    "Orchidaceae_12",
    "Orchidaceae_13",
    "Orchidaceae_14",
    "Orchidaceae_15",
    "Orchidaceae_16",
    "Orchidaceae_17",
    "Orchidaceae_18",
    "Orchidaceae_19",
    "Orchidaceae_1",
    "Orchidaceae_20",
    "Orchidaceae_21",
    "Orchidaceae_22",
    "Orchidaceae_23",
    "Orchidaceae_24",
    "Orchidaceae_25",
    "Orchidaceae_2",
    "Orchidaceae_3",
    "Orchidaceae_4",
    "Orchidaceae_5",
    "Orchidaceae_6",
    "Orchidaceae_7",
    "Orchidaceae_8",
    "Orchidaceae_9",
    "Orobanchaceae_Phrymaceae_Mazaceae_Paulowniaceae_1",
    "Orobanchaceae_Phrymaceae_Mazaceae_Paulowniaceae_2",
    "Orobanchaceae_Phrymaceae_Mazaceae_Paulowniaceae_3",
    "Orobanchaceae_Phrymaceae_Mazaceae_Paulowniaceae_4",
    "Orobanchaceae_Phrymaceae_Mazaceae_Paulowniaceae_5",
    "Papaveraceae_1",
    "Papaveraceae_2",
    "Phyllanthaceae_Picodendraceae_Linaceae_Ixonanthaceae_1",
    "Phyllanthaceae_Picodendraceae_Linaceae_Ixonanthaceae_2",
    "Phyllanthaceae_Picodendraceae_Linaceae_Ixonanthaceae_3",
    "Phyllanthaceae_Picodendraceae_Linaceae_Ixonanthaceae_4",
    "Phyllanthaceae_Picodendraceae_Linaceae_Ixonanthaceae_5",
    "Phyllanthaceae_Picodendraceae_Linaceae_Ixonanthaceae_6",
    "Plantaginaceae_1",
    "Plantaginaceae_2",
    "Plantaginaceae_3",
    "Plumbaginaceae_Polygonaceae_Frankeniaceae_Tamaricaceae_1",
    "Plumbaginaceae_Polygonaceae_Frankeniaceae_Tamaricaceae_2",
    "Plumbaginaceae_Polygonaceae_Frankeniaceae_Tamaricaceae_3",
    "Plumbaginaceae_Polygonaceae_Frankeniaceae_Tamaricaceae_4",
    "Poaceae_10",
    "Poaceae_11",
    "Poaceae_12",
    "Poaceae_13",
    "Poaceae_14",
    "Poaceae_15",
    "Poaceae_16",
    "Poaceae_17",
    "Poaceae_1",
    "Poaceae_2",
    "Poaceae_3",
    "Poaceae_4",
    "Poaceae_5",
    "Poaceae_6",
    "Poaceae_7",
    "Poaceae_8",
    "Poaceae_9",
    "Ranunculaceae_1",
    "Ranunculaceae_2",
    "Ranunculaceae_3",
    "Rosaceae_1",
    "Rosaceae_2",
    "Rosaceae_3",
    "Rosaceae_4",
    "Rosaceae_5",
    "Rosaceae_6",
    "Rubiaceae_1",
    "Rubiaceae_2",
    "Rubiaceae_3",
    "Scrophulariaceae_1",
    "Zingiberaceae_1",
    "Zingiberaceae_2",
    "Aizoaceae_Phytolaccaceae_Barbeuiaceae_Lophiocarpaceae_Gisekiaceae_Sarcobataceae",
    "Calyceraceae",
    "Capparaceae",
    "Chrysobalanaceae_Malpighiaceae_Caryocaraceae_Balanopaceae_Elatinaceae_Centroplacaceae_Dichapetalaceae_Putranjivaceae_Euphroniaceae_Lophopyxidaceae_Trigoniaceae",
    "Cleomaceae",
    "Crassulaceae_Aphanopetalaceae_Halograceae_Penthoraceae_Tetracarpaeaceae",
    "Dipterocarpaceae_Bixaceae_Cistaceae_Sarcoleanaceae_Muntingiaceae_Sphaerosepalaceae",
    "Simaroubaceae",
    "Alzateaceae_Crypteroniaceae_Penaeaceae",
    "Araliaceae", 
    "Asphodelaceae",
    "Balsaminaceae_Marcgraviaceae_Tetrameristaceae",
    "Berberidaceae", 
    "Bignoniaceae",
    "Cactaceae_Molluginaceae_Didiereaceae_Anacompserotaceae_Basellaceae_Montiaceae_Halophytaceae_Portulacaceae_Talinaceae",
    "Campanulaceae_Rousseaceae",
    "Cannabaceae",
    "Cercidiphyllaceae_Hamamelidaceae_Daphniphyllaceae_Altingiaceae_Paeoniaceae", 
    "Combretaceae", 
    "Droseraceae_Ancistrocladaceae_Drosophyllaceae_Nepenthaceae_Dioncophyllaceae", 
    "Ebenaceae", 
    "Gentianaceae", 
    "Goodeniaceae", 
    "Juncaceae", 
    "Loganiaceae_Gelsemiaceae", 
    "Lythraceae_Onagraceae", 
    "Malvaceae", 
    "Meliaceae", 
    "Menispermaceae", 
    "Menyanthaceae", 
    "Monimiaceae", 
    "Moraceae", 
    "Passifloraceae", 
    "Pentaphylacaceae_Sladeniaceae", 
    "Pittosporaceae", 
    "Polemoniaceae_Lecythidaceae_Fouquieriaceae", 
    "Polygalaceae_Surianaceae", 
    "Primulaceae", 
    "Resedaceae", 
    "Restionaceae", 
    "Rhamnaceae_Barbeyaceae_Dirachmaceae_Elaeagnaceae", 
    "Rutaceae", 
    "Salicaceae_Lacistemataceae", 
    "Sapindaceae", 
    "Sapotaceae", 
    "Saxifragaceae_Iteaceae_Grossulariaceae", 
    "Styracaceae_Diapensiaceae_Symplocaceae", 
    "Theaceae", 
    "Thymelaeaceae", 
    "Typhaceae", 
    "Ulmaceae", 
    "Urticaceae", 
    "Verbenaceae_Schlegeliaceae_Lentibulariaceae_Thomandersiaceae", 
    "Violaceae_Goupiaceae", 
    "Xyridaceae_Eriocaulaceae", 
    "Buxales",
    "Canellales",
    "Celastrales",
    "Commelinales",
    "Cornales",
    "Cupressales",
    "Escalloniales",
    "Fagales",
    "Geraniales",
    "Gnetales",
    "Huerteales",
    "Icacinales",
    "Metteniusales",
    "Nymphaeales",
    "Oxalidales",
    "Pandanales",
    "Pinales",
    "Piperales",
    "Zygophyllales")

# Select the elements from file list which are in the orders_fams_subfams_list
file_list_subset <- file_list[grep(paste(orders_fams_subfams_list, collapse = "|"), file_list)]
file_list_subset

# Now we load all these trees into a list of phylogenies
phylogenies <- lapply(file_list_subset, read.tree)

# Now can we find the number of species across all these trees
Nr_species_esse <- sum(sapply(phylogenies, function(x) length(x$tip.label)))
Nr_species_esse # 61291
67379-61291

# What is the Mean min and max number of species per tree?
mean_species_per_tree <- mean(sapply(phylogenies, function(x) length(x$tip.label)))
mean_species_per_tree # 182

# What is the Mean min and max number of species per tree?
min_species_per_tree <- min(sapply(phylogenies, function(x) length(x$tip.label)))
min_species_per_tree # 4

# Max species per tree
max_species_per_tree <- max(sapply(phylogenies, function(x) length(x$tip.label)))
max_species_per_tree # 981

# How many phylogenies have less than 100 species?
Nr_phylogenies_less_than_100 <- sum(sapply(phylogenies, function(x) length(x$tip.label)) < 100)
Nr_phylogenies_less_than_100 # 138

# How many phylogenies have less than 10 species?
Nr_phylogenies_less_than_10 <- sum(sapply(phylogenies, function(x) length(x$tip.label)) < 10)
Nr_phylogenies_less_than_10 # 2
# which phylogenies have less than 10 species?
phylogenies[sapply(phylogenies, function(x) length(x$tip.label)) < 10]
