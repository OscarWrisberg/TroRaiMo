################################################################################################################################################
######################################################-- Loading packages --####################################################################
################################################################################################################################################
# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("data.table", "ape", "phytools", "geiger", "castor", "MonoPhy", "terra", "dplyr", "ggplot2", "sp","sf")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))



################################################################################################################################################
#######################################################-- Local testing --######################################################################
################################################################################################################################################
# srun --account Trf_models --mem 100g --pty bash
#setwd("/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/orders") # local you need around 60 gigs of ram
setwd("/home/owrisberg/Trf_models/workflow/02_adding_orders/pruning/subset_of_orders") # srun
input_file_tree <- "sub_phylo_Apiaceae_2.tre"
output <- "Test_Apiaceae_2.txt"
input_file_wcvp <- "/home/owrisberg/Trf_models/workflow/02_adding_orders/wcvp_names_apg_aligned.rds" #srun
path_out <- "/home/owrisberg/Trf_models//workflow/03_distribution_data/" #srun
order_in_question <- as.character("Apiaceae_2")
apg  <- "../../../../TroRaiMo/apgweb_parsed.csv"
renamed_occurence_file <- "/home/owrisberg/Trf_models//workflow/01_distribution_data/06_Renamed/gbif_renamed.rds" #srun
koppen_biome_file <- "../../../../TroRaiMo/koppen_geiger_0p01.tif" #srun
percentages_for_present <- 0.33

################################################################################################################################################
##############################################-- Handling Command Line arguments --#############################################################
################################################################################################################################################

# Command line arguments
input_file_tree <- commandArgs(trailingOnly = TRUE)[1]
output <- commandArgs(trailingOnly = TRUE)[2]
input_file_wcvp <- commandArgs(trailingOnly = TRUE)[3]
path_out <- commandArgs(trailingOnly = TRUE)[4]
order_in_question <- commandArgs(trailingOnly = TRUE)[5]
apg <- commandArgs(trailingOnly = TRUE)[6]
renamed_occurence_file <- commandArgs(trailingOnly = TRUE)[7]
koppen_biome_file <- commandArgs(trailingOnly = TRUE)[8]
percentages_for_present <- commandArgs(trailingOnly = TRUE)[9]

# Print the command line arguments
cat("The input file for the tree is ", input_file_tree, "\n")
cat("The output file is ", output, "\n")
cat("The input file for the wcvp dataset is ", input_file_wcvp, "\n")
cat("The path to the output file is ", path_out, "\n")
cat("The order of the tree is ", order_in_question, "\n")
cat("The current working directory is ", getwd(), "\n")

################################################################################################################################################
######################################################-- Loading files --#######################################################################
################################################################################################################################################

# Load the wcvp dataset
wcvp <- readRDS(input_file_wcvp)

# Load the tree
tree <- read.tree(input_file_tree)

# Load apg
apg <- fread(apg)

#Load renamed_occurence
renamed_occurence <- readRDS(renamed_occurence_file)

# Filter the wcvp dataset to include only rows where taxon_status == "Accepted"
wcvp_accepted <- subset(wcvp, taxon_status == "Accepted")
wcvp_accepted_species <- subset(wcvp_accepted, taxon_rank == "Species")
wcvp_accepted_species <- subset(wcvp_accepted_species, species_hybrid == "")

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

################################################################################################################################################
#############################################-- Finding Environmental data GBIF --##############################################################
################################################################################################################################################

# Subset the renamed_occurence dataset to only include the genera in the tree
cat("Subsetting the renamed_occurence dataset to only include the genera in the tree \n\n")

# Finding the number of species in the tree based on the 
# Splitting all the tip labels into genus and species
tree_genera <- sapply(strsplit(tree$tip.label, " "), "[", 1)
unique_genera_in_tree <- unique(tree_genera)

# Subsetting the renamed occurrence dataset to only include the genera in the tree
renamed_occurence_subset <- renamed_occurence[which(renamed_occurence$genus %in% unique_genera_in_tree),]
cat("There are ",length(unique(renamed_occurence_subset$species)), " species in the renamed_occurence_subset \n\n") 

#Remove occurences with NA in decimalLatitude or decimalLongitude
#cat("Removing occurences with NA in decimalLatitude or decimalLongitude \n\n")
renamed_occurence_subset <- renamed_occurence_subset[which(!is.na(renamed_occurence_subset$decimalLatitude) & !is.na(renamed_occurence_subset$decimalLongitude)),]


# There is an option to add additional cleaning of the data here, but for now we will continue with the data as it is.

# How many of the species in the tree are found in the renamed_occurence_subset
cat("Out of ",length(tree$tip.label),"tips in the tree there are ",length(tree$tip.label[which(tree$tip.label %in% renamed_occurence_subset$wcvp_taxon_name)]),"of the species in the tree are found in the renamed_occurence_subset \n\n")

# Loading the Koppen biomes data
cat("Loading the Koppen biomes data \n\n")
koppen_biome_map <- rast(koppen_biome_file)

##################################################################################################################################################################
############################################--- Find the estimated biome based on occurrences ---#################################################################
##################################################################################################################################################################

# Find all the genera which are found in the tree.
unique_genera_in_tree <- unique(wcvp_accepted_species$genus[which(wcvp_accepted_species$taxon_name %in% tree$tip.label)])
unique_genera_in_tree

# Assign Koppen biomes to Tropical rainforest or not.
cat("Adding a new column to the csv file which is Tropical rainforest or not \n\n")
koppen_biome_map_tropical_rainforest <- ifel(koppen_biome_map == 1 | koppen_biome_map == 2, 1, 0)

# Extracting the biome for each of the occurrences
renamed_occurence_subset$in_tropical_rainforest <- terra::extract(koppen_biome_map_tropical_rainforest, renamed_occurence_subset[,c("decimalLongitude", "decimalLatitude")])[,2]

# Remove occurences with NA in in_tropical_rainforest
cat("Removing ", length(which(is.na(terra::extract(koppen_biome_map_tropical_rainforest, renamed_occurence_subset[,c("decimalLongitude", "decimalLatitude")])[,2]))) ," occurences with NA in in_tropical_rainforest \n\n")
renamed_occurence_subset <- renamed_occurence_subset[which(!is.na(renamed_occurence_subset$in_tropical_rainforest)),]

# Summarize the number of occurrences which are inside and outside the tropical rainforest
result_summary <- aggregate(renamed_occurence_subset$in_tropical_rainforest, by = list(renamed_occurence_subset$wcvp_taxon_name), FUN = function(x) c(sum(x == 1), sum(x == 0)))

# Splitting the second column, which is a column containing a dataframe with 2 columns, into two individual columns
result_summary <- cbind(result_summary[,1] , as.data.frame(result_summary[,2]))

# Changing the column names
result_summary <- setNames(result_summary, c("wcvp_taxon_name", "occurrences_trf", "occurrences_non_trf"))

# Calculate the proportion of occurrences inside the tropical rainforest biome
result_summary$proportion_in_tropical_rainforest <- result_summary$occurrences_trf/(result_summary$occurrences_trf + result_summary$occurrences_non_trf)
result_summary$proportion_outside_tropical_rainforest <- result_summary$occurrences_non_trf/(result_summary$occurrences_trf + result_summary$occurrences_non_trf)

# Now we merge the result_summary with the wcvp_accepted_species to get the climate description for each species
result_summary <- merge(result_summary, wcvp_accepted_species[,c("taxon_name", "climate_description")], by.x = "wcvp_taxon_name", by.y = "taxon_name", all.x = TRUE)
head(result_summary)

#How many of the species in result_summary have the climate description wet_tropical.
cat("Out of ", nrow(result_summary), " species in result_summary ", length(which(result_summary$climate_description == "wet tropical")), " have the climate description wet tropical \n\n")

# We now calculate if species are Trf, non-trf or widespread based on the percentages for present
result_summary$trf <- ifelse(result_summary$proportion_in_tropical_rainforest >= percentages_for_present, 1, 0)
result_summary$non_trf <- ifelse(result_summary$proportion_outside_tropical_rainforest >= percentages_for_present, 1, 0)

# How many percent of the wet tropical species are Trf
Wet_tropical_trf <- length(which(result_summary$trf == 1 & result_summary$climate_description == "wet tropical" & result_summary$non_trf == 0))
Wet_tropical_trf_percentage <- Wet_tropical_trf/length(which(result_summary$climate_description == "wet tropical"))
if(is.nan(Wet_tropical_trf_percentage)){
  Wet_tropical_trf_percentage <- 0
}


# How many percent of the wet tropical species are widespread
Wet_tropical_widespread <- length(which(result_summary$trf == 1 & result_summary$climate_description == "wet tropical" & result_summary$non_trf == 1))
Wet_tropical_widespread_percentage <- Wet_tropical_widespread/length(which(result_summary$climate_description == "wet tropical"))
if(is.nan(Wet_tropical_widespread_percentage)){
  Wet_tropical_widespread_percentage <- 0
}

# How many percent of the wet tropical species are non-trf
Wet_tropical_non_trf <- length(which(result_summary$trf == 0 & result_summary$climate_description == "wet tropical" & result_summary$non_trf == 1))
Wet_tropical_non_trf_percentage <- Wet_tropical_non_trf/length(which(result_summary$climate_description == "wet tropical"))
if(is.nan(Wet_tropical_non_trf_percentage)){
  Wet_tropical_non_trf_percentage <- 0
}

# How many percent of the non-wet tropical species are Trf
non_wet_tropical_trf <- length(which(result_summary$trf == 1 & result_summary$climate_description != "wet tropical" & result_summary$non_trf == 0))
non_wet_tropical_trf_percentage <- non_wet_tropical_trf/length(which(result_summary$climate_description != "wet tropical"))
if(is.nan(non_wet_tropical_trf_percentage)){
  non_wet_tropical_trf_percentage <- 0
}

# How many percent of the non-wet tropical species are widespread
non_wet_tropical_widespread <- length(which(result_summary$trf == 1 & result_summary$climate_description != "wet tropical" & result_summary$non_trf == 1))
non_wet_tropical_widespread_percentage <- non_wet_tropical_widespread/length(which(result_summary$climate_description != "wet tropical"))
if(is.nan(non_wet_tropical_widespread_percentage)){
  non_wet_tropical_widespread_percentage <- 0
}

# How many percent of the non-wet tropical species are non-trf
non_wet_tropical_non_trf <- length(which(result_summary$trf == 0 & result_summary$climate_description != "wet tropical" & result_summary$non_trf == 1))
non_wet_tropical_non_trf_percentage <- non_wet_tropical_non_trf/length(which(result_summary$climate_description != "wet tropical"))
if(is.nan(non_wet_tropical_non_trf_percentage)){
  non_wet_tropical_non_trf_percentage <- 0
}

#I therefore need find all the species which are not in the occurrences but are in wcvp
wcvp_accepted_species_subset <- wcvp_accepted_species[which(wcvp_accepted_species$genus %in% unique_genera_in_tree),]

##################################################################################################################################################################
##############################################--- Find the number of species per biome for the clade ---##########################################################
##################################################################################################################################################################

# First we find the biome of all the species which are in occurrences.
wcvp_clade <- wcvp_accepted_species_subset[which(wcvp_accepted_species_subset$genus %in% unique_genera_in_tree),] # getting all the accepted species from wcvp which are in the genera in the tree

# The biome of all the species which are in occurrences
result_summary_total_sp <- result_summary

# Finding the species which are not in the occurrences
Sp_not_not_in_occurrences <- setdiff(wcvp_accepted_species_subset$taxon_name, result_summary_total_sp$wcvp_taxon_name)

# Then we find their biome using wcvp
wcvp_clade_not_in_occurrences <- wcvp_accepted_species_subset[which(wcvp_accepted_species_subset$taxon_name %in% Sp_not_not_in_occurrences),]

# Find the number of species which are wet tropical and non-wet tropical
nr_wet_tropical <- length(which(wcvp_clade_not_in_occurrences$climate_description == "wet tropical"))
nr_non_wet_tropical <- length(which(wcvp_clade_not_in_occurrences$climate_description != "wet tropical"))

# We now multiply these numbers with our estimates on Trf, non_trf and widespread
nr_wet_tropical_trf <- nr_wet_tropical * Wet_tropical_trf_percentage
nr_wet_tropical_non_trf <- nr_wet_tropical * Wet_tropical_non_trf_percentage
nr_wet_tropical_widespread <- nr_wet_tropical * Wet_tropical_widespread_percentage

nr_non_wet_tropical_trf <- nr_non_wet_tropical * non_wet_tropical_trf_percentage
nr_non_wet_tropical_non_trf <- nr_non_wet_tropical * non_wet_tropical_non_trf_percentage
nr_non_wet_tropical_widespread <- nr_non_wet_tropical * non_wet_tropical_widespread_percentage

# Find the number of Trf, non_trf and widespread species in result_summary_total_sp
nr_trf <- length(which(result_summary_total_sp$trf == 1))
nr_non_trf <- length(which(result_summary_total_sp$non_trf == 1))
nr_widespread <- length(which(result_summary_total_sp$trf == 1 & result_summary_total_sp$non_trf == 1))

# Adding these numbers together
nr_trf <- nr_trf + nr_wet_tropical_trf + nr_non_wet_tropical_trf
nr_non_trf <- nr_non_trf + nr_wet_tropical_non_trf + nr_non_wet_tropical_non_trf
nr_widespread <- nr_widespread + nr_wet_tropical_widespread + nr_non_wet_tropical_widespread


# So now we have the states for all the tips in the tree
# Then we find the states for all the species not in the tree.
sp_not_in_tree <- setdiff(wcvp_accepted_species_subset$taxon_name, tree$tip.label)


##################################################################################################################################################################
##########################################################--- Find the estimated biome based on occurrences ---###################################################
##################################################################################################################################################################

if (all(sp_not_in_tree %in% renamed_occurence_subset$wcvp_taxon_name)) { # If all the tip labels are found in the occurences use only them to create the presence absence matrix

  # Subset the occurrences to only include the species that are missing from the tree
  renamed_occurence_subset_missing_sp <- renamed_occurence[which(renamed_occurence$wcvp_taxon_name %in% sp_not_in_tree),]

  cat("Adding a new column to the csv file which is Tropical rainforest or not \n\n")
  koppen_biome_map_tropical_rainforest <- ifel(koppen_biome_map == 1 | koppen_biome_map == 2, 1, 0)

  # Extracting the biome for each of the occurrences
  renamed_occurence_subset$in_tropical_rainforest <- terra::extract(koppen_biome_map_tropical_rainforest, renamed_occurence_subset[,c("decimalLongitude", "decimalLatitude")])[,2]

  # Remove occurences with NA in in_tropical_rainforest
  cat("Removing ", length(which(is.na(terra::extract(koppen_biome_map_tropical_rainforest, renamed_occurence_subset[,c("decimalLongitude", "decimalLatitude")])[,2]))) ," occurences with NA in in_tropical_rainforest \n\n")
  renamed_occurence_subset <- renamed_occurence_subset[which(!is.na(renamed_occurence_subset$in_tropical_rainforest)),]

  # Summarize the number of occurrences which are inside and outside the tropical rainforest
  missing_sp_result_summary <- aggregate(renamed_occurence_subset$in_tropical_rainforest, by = list(renamed_occurence_subset$wcvp_taxon_name), FUN = function(x) c(sum(x == 1), sum(x == 0)))

  # Calculate the proportion of occurrences inside the tropical rainforest biome
  missing_sp_result_summary$proportion_in_tropical_rainforest <- missing_sp_result_summary$x[,1]/(missing_sp_result_summary$x[,1] + missing_sp_result_summary$x[,2])
  missing_sp_result_summary$proportion_outside_tropical_rainforest <- missing_sp_result_summary$x[,2]/(missing_sp_result_summary$x[,1] + missing_sp_result_summary$x[,2])

  # Use this summary to determine how many of the missing species are Trf, non-trf or widespread
  missing_sp_trf <- length(which(missing_sp_result_summary$proportion_in_tropical_rainforest >= percentages_for_present))
  missing_sp_non_trf <- length(which(missing_sp_result_summary$proportion_outside_tropical_rainforest >= percentages_for_present))
  missing_sp_widespread <- length(which(missing_sp_result_summary$proportion_in_tropical_rainforest < percentages_for_present & missing_sp_result_summary$proportion_outside_tropical_rainforest < percentages_for_present))
  missing_sp_biomes <- c(missing_sp_trf, missing_sp_non_trf, missing_sp_widespread)

  # Find the number of 

  #return(missing_sp_biomes)
  missing_sp_biomes

##################################################################################################################################################################
######################################################--- Find the species which are not in occurrences  ---######################################################
##################################################################################################################################################################

} else { # If all the species are not found in the occurrences, see if the species missing from the occurrences are found in the WCVP
  renamed_occurence_subset <- renamed_occurence[which(renamed_occurence$wcvp_taxon_name %in% sp_not_in_tree),] # Subsetting the renamed_occurence dataset to only include the species in the tree
  cat("Not all the species not in the tree are found in the renamed_occurence_subset \n\n")
  cat("The following species are not found in the renamed_occurence_subset: \n")
  cat(setdiff(sp_not_in_tree, renamed_occurence_subset$wcvp_taxon_name), "\n\n")
  occurrence_missing_sp <- setdiff(sp_not_in_tree, renamed_occurence_subset$wcvp_taxon_name) # Finding the species which are missing from the occurrences.


  cat("Can we find a biome for the missing species in the wcvp dataset? \n\n")
  # subsetting the wcvp dataset to only include the species in the tree
  wcvp_subset <- wcvp[which(wcvp$taxon_name %in% sp_not_in_tree),]

  # Subsetting wcvp to only include accepted species and species with a climate description
  wcvp_subset <- wcvp_subset[which(wcvp_subset$taxon_status == "Accepted" & wcvp_subset$climate_description != "" & wcvp_subset$taxon_rank == "Species" & wcvp_subset$species_hybrid == ""),]
  cat("Are all the missing sp found in WCVP and do they have a Climate description: ", all(occurrence_missing_sp %in% wcvp_subset$taxon_name), "\n")

##################################################################################################################################################################
###############################################--- Find the estimated biome for sp not in occurrences based on wcvp ---###########################################
##################################################################################################################################################################

  if (all(occurrence_missing_sp %in% wcvp_subset$taxon_name)) { # If all the missing species are found in the wcvp use the wcvp climate column to add the missing species to the presence absence matrix

    # use the koppen biome maps for the species represented in the renamed_occurence_subset
    # Add a new column to the csv file which is Tropical rainforest or not.
    cat("Adding a new column to the csv file which is Tropical rainforest or not \n\n")
    koppen_biome_map_tropical_rainforest <- ifel(koppen_biome_map == 1 | koppen_biome_map == 2, 1, 0)

    # Extracting the biome for each of the occurrences
    renamed_occurence_subset$in_tropical_rainforest <- terra::extract(koppen_biome_map_tropical_rainforest, renamed_occurence_subset[,c("decimalLongitude", "decimalLatitude")])[,2]

    # Remove occurences with NA in in_tropical_rainforest
    cat("Removing ", length(which(is.na(terra::extract(koppen_biome_map_tropical_rainforest, renamed_occurence_subset[,c("decimalLongitude", "decimalLatitude")])[,2]))) ," occurences with NA in in_tropical_rainforest \n\n")
    renamed_occurence_subset <- renamed_occurence_subset[which(!is.na(renamed_occurence_subset$in_tropical_rainforest)),]

    # Summarize the number of occurrences which are inside and outside the tropical rainforest
    result_summary_missing_sp <- aggregate(renamed_occurence_subset$in_tropical_rainforest, by = list(renamed_occurence_subset$wcvp_taxon_name), FUN = function(x) c(sum(x == 1), sum(x == 0)))

    # Splitting the second column, which is a column containing a dataframe with 2 columns, into two individual columns
    result_summary_missing_sp <- cbind(result_summary_missing_sp[,1] , as.data.frame(result_summary_missing_sp[,2]))
    
    # Changing the column names
    result_summary_missing_sp <- setNames(result_summary_missing_sp, c("wcvp_taxon_name", "occurrences_trf", "occurrences_non_trf"))

    # Calculate the proportion of occurrences inside the tropical rainforest biome
    result_summary_missing_sp$proportion_in_tropical_rainforest <- result_summary_missing_sp$occurrences_trf/(result_summary_missing_sp$occurrences_trf + result_summary_missing_sp$occurrences_non_trf)
    result_summary_missing_sp$proportion_outside_tropical_rainforest <- result_summary_missing_sp$occurrences_non_trf/(result_summary_missing_sp$occurrences_trf + result_summary_missing_sp$occurrences_non_trf)

    # Now we sum up how many of the missing species are Trf, non-trf or widespread
    missing_sp_trf <- length(which(result_summary_missing_sp$proportion_in_tropical_rainforest >= percentages_for_present))
    missing_sp_non_trf <- length(which(result_summary_missing_sp$proportion_outside_tropical_rainforest >= percentages_for_present))
    missing_sp_widespread <- length(which(result_summary_missing_sp$proportion_in_tropical_rainforest < percentages_for_present & result_summary_missing_sp$proportion_outside_tropical_rainforest < percentages_for_present))

    # adding them together to a vector
    missing_sp_biomes <- c(missing_sp_trf, missing_sp_non_trf, missing_sp_widespread)

    # We then need to add the counts from the 
    # Use the wcvp_subset to find the biome for the rest of the species
    cat("Using the wcvp_subset to find the biome for the rest of the species \n\n")

    # I will then find their climate description in the wcvp_subset
    missing_sp_climate <- wcvp_subset[which(wcvp_subset$taxon_name %in% occurrence_missing_sp),]

    # We will now find the number of the missing species which are wet tropical and non-wet tropical
    nr_wet_tropical <- length(which(missing_sp_climate$climate_description == "wet tropical"))
    nr_non_wet_tropical <- length(which(missing_sp_climate$climate_description != "wet tropical"))

    # We now multiply these numbers with the percentages to find out how many of the missing species are Trf, non-trf or widespread
    missing_sp_trf_wcvp <- nr_wet_tropical * Wet_tropical_trf_percentage + nr_non_wet_tropical * non_wet_tropical_trf_percentage
    missing_sp_non_trf_wcvp <- nr_wet_tropical * Wet_tropical_non_trf_percentage + nr_non_wet_tropical * non_wet_tropical_non_trf_percentage
    missing_sp_widespread_wcvp <- nr_wet_tropical * Wet_tropical_widespread_percentage + nr_non_wet_tropical * non_wet_tropical_widespread_percentage

    # adding them to the missing_sp_biomes vector
    missing_sp_biomes[1] <- missing_sp_biomes[1] + missing_sp_trf_wcvp
    missing_sp_biomes[2] <- missing_sp_biomes[2] + missing_sp_non_trf_wcvp
    missing_sp_biomes[3] <- missing_sp_biomes[3] + missing_sp_widespread_wcvp

    #return(missing_sp_biomes)
    missing_sp_biomes

##################################################################################################################################################################
#########################################################--- Removing species not in occurrences or wcvp ---######################################################
##################################################################################################################################################################

  } else {
    missing_sp <- occurrence_missing_sp[which(!occurrence_missing_sp %in% wcvp_subset$taxon_name)]
    cat("Not all the missing species are found in the wcvp dataset \n")
    cat("The following species are found in the wcvp dataset \n")
    cat(missing_sp, "\n\n")
    cat("The following species are not found in the wcvp dataset \n")
    cat(setdiff(missing_sp, wcvp_subset$taxon_name), "\n\n")
    wcvp_missing_sp <- setdiff(missing_sp, wcvp_subset$taxon_name)

    missing_sp_file_name <- paste0("Missing_sp_",order_in_question,".txt")

    cat("In order for the pipeline to progess were removing the missing species from the tree, but they can be found in a file called ", missing_sp_file_name, "\n")
    write.table(wcvp_missing_sp, file = paste0(path_out, "Missing_sp_for_biome_sampling_fraction",order_in_question,".txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

    # dropping the missing species from the sp_not_in_tree.
    sp_not_in_tree <- setdiff(sp_not_in_tree, wcvp_missing_sp)

####################################################################################################################################################################
################################################### --- If all species are in occurrences after the removal ---#####################################################
####################################################################################################################################################################
  if (all(sp_not_in_tree %in% renamed_occurence_subset$wcvp_taxon_name)){

    # Subset the occurrences to only include the species that are missing from the tree
    renamed_occurence_subset_missing_sp <- renamed_occurence[which(renamed_occurence$wcvp_taxon_name %in% sp_not_in_tree),]

    cat("Adding a new column to the csv file which is Tropical rainforest or not \n\n")
    koppen_biome_map_tropical_rainforest <- ifel(koppen_biome_map == 1 | koppen_biome_map == 2, 1, 0)

    # Extracting the biome for each of the occurrences
    renamed_occurence_subset$in_tropical_rainforest <- terra::extract(koppen_biome_map_tropical_rainforest, renamed_occurence_subset[,c("decimalLongitude", "decimalLatitude")])[,2]

    # Remove occurences with NA in in_tropical_rainforest
    cat("Removing ", length(which(is.na(terra::extract(koppen_biome_map_tropical_rainforest, renamed_occurence_subset[,c("decimalLongitude", "decimalLatitude")])[,2]))) ," occurences with NA in in_tropical_rainforest \n\n")
    renamed_occurence_subset <- renamed_occurence_subset[which(!is.na(renamed_occurence_subset$in_tropical_rainforest)),]

    # Summarize the number of occurrences which are inside and outside the tropical rainforest
    missing_sp_result_summary <- aggregate(renamed_occurence_subset$in_tropical_rainforest, by = list(renamed_occurence_subset$wcvp_taxon_name), FUN = function(x) c(sum(x == 1), sum(x == 0)))

    # Calculate the proportion of occurrences inside the tropical rainforest biome
    missing_sp_result_summary$proportion_in_tropical_rainforest <- missing_sp_result_summary$x[,1]/(missing_sp_result_summary$x[,1] + missing_sp_result_summary$x[,2])
    missing_sp_result_summary$proportion_outside_tropical_rainforest <- missing_sp_result_summary$x[,2]/(missing_sp_result_summary$x[,1] + missing_sp_result_summary$x[,2])

    # Use this summary to determine how many of the missing species are Trf, non-trf or widespread
    missing_sp_trf <- length(which(missing_sp_result_summary$proportion_in_tropical_rainforest >= percentages_for_present))
    missing_sp_non_trf <- length(which(missing_sp_result_summary$proportion_outside_tropical_rainforest >= percentages_for_present))
    missing_sp_widespread <- length(which(missing_sp_result_summary$proportion_in_tropical_rainforest < percentages_for_present & missing_sp_result_summary$proportion_outside_tropical_rainforest < percentages_for_present))
    missing_sp_biomes <- c(missing_sp_trf, missing_sp_non_trf, missing_sp_widespread)

    #return(missing_sp_biomes)
    missing_sp_biomes

####################################################################################################################################################################
############################################## --- Find the estimated biome for sp not in occurrences based on wcvp ---#############################################
####################################################################################################################################################################

  } else {
        # use the koppen biome maps for the species represented in the renamed_occurence_subset
    # Add a new column to the csv file which is Tropical rainforest or not.
    cat("Adding a new column to the csv file which is Tropical rainforest or not \n\n")
    koppen_biome_map_tropical_rainforest <- ifel(koppen_biome_map == 1 | koppen_biome_map == 2, 1, 0)

    # Extracting the biome for each of the occurrences
    renamed_occurence_subset$in_tropical_rainforest <- terra::extract(koppen_biome_map_tropical_rainforest, renamed_occurence_subset[,c("decimalLongitude", "decimalLatitude")])[,2]

    # Remove occurences with NA in in_tropical_rainforest
    cat("Removing ", length(which(is.na(terra::extract(koppen_biome_map_tropical_rainforest, renamed_occurence_subset[,c("decimalLongitude", "decimalLatitude")])[,2]))) ," occurences with NA in in_tropical_rainforest \n\n")
    renamed_occurence_subset <- renamed_occurence_subset[which(!is.na(renamed_occurence_subset$in_tropical_rainforest)),]

    # Summarize the number of occurrences which are inside and outside the tropical rainforest
    result_summary_missing_sp <- aggregate(renamed_occurence_subset$in_tropical_rainforest, by = list(renamed_occurence_subset$wcvp_taxon_name), FUN = function(x) c(sum(x == 1), sum(x == 0)))

    # Splitting the second column, which is a column containing a dataframe with 2 columns, into two individual columns
    result_summary_missing_sp <- cbind(result_summary_missing_sp[,1] , as.data.frame(result_summary_missing_sp[,2]))
    
    # Changing the column names
    result_summary_missing_sp <- setNames(result_summary_missing_sp, c("wcvp_taxon_name", "occurrences_trf", "occurrences_non_trf"))

    # Calculate the proportion of occurrences inside the tropical rainforest biome
    result_summary_missing_sp$proportion_in_tropical_rainforest <- result_summary_missing_sp$occurrences_trf/(result_summary_missing_sp$occurrences_trf + result_summary_missing_sp$occurrences_non_trf)
    result_summary_missing_sp$proportion_outside_tropical_rainforest <- result_summary_missing_sp$occurrences_non_trf/(result_summary_missing_sp$occurrences_trf + result_summary_missing_sp$occurrences_non_trf)

    # Now we sum up how many of the missing species are Trf, non-trf or widespread
    missing_sp_trf <- length(which(result_summary_missing_sp$proportion_in_tropical_rainforest >= percentages_for_present))
    missing_sp_non_trf <- length(which(result_summary_missing_sp$proportion_outside_tropical_rainforest >= percentages_for_present))
    missing_sp_widespread <- length(which(result_summary_missing_sp$proportion_in_tropical_rainforest < percentages_for_present & result_summary_missing_sp$proportion_outside_tropical_rainforest < percentages_for_present))

    # adding them together to a vector
    missing_sp_biomes <- c(missing_sp_trf, missing_sp_non_trf, missing_sp_widespread)

    # We then need to add the counts from the 
    # Use the wcvp_subset to find the biome for the rest of the species
    cat("Using the wcvp_subset to find the biome for the rest of the species \n\n")

    # I will then find their climate description in the wcvp_subset
    missing_sp_climate <- wcvp_subset[which(wcvp_subset$taxon_name %in% occurrence_missing_sp),]

    # We will now find the number of the missing species which are wet tropical and non-wet tropical
    nr_wet_tropical <- length(which(missing_sp_climate$climate_description == "wet tropical"))
    nr_non_wet_tropical <- length(which(missing_sp_climate$climate_description != "wet tropical"))

    # We now multiply these numbers with the percentages to find out how many of the missing species are Trf, non-trf or widespread
    missing_sp_trf_wcvp <- nr_wet_tropical * Wet_tropical_trf_percentage + nr_non_wet_tropical * non_wet_tropical_trf_percentage
    missing_sp_non_trf_wcvp <- nr_wet_tropical * Wet_tropical_non_trf_percentage + nr_non_wet_tropical * non_wet_tropical_non_trf_percentage
    missing_sp_widespread_wcvp <- nr_wet_tropical * Wet_tropical_widespread_percentage + nr_non_wet_tropical * non_wet_tropical_widespread_percentage

    # adding them to the missing_sp_biomes vector
    missing_sp_biomes[1] <- missing_sp_biomes[1] + missing_sp_trf_wcvp
    missing_sp_biomes[2] <- missing_sp_biomes[2] + missing_sp_non_trf_wcvp
    missing_sp_biomes[3] <- missing_sp_biomes[3] + missing_sp_widespread_wcvp

    #return(missing_sp_biomes)
    missing_sp_biomes

    }
  }
}

# Find the total number of accepted species in the genera which are in the tree wcvp
total_number_accepted_species <- length(wcvp_accepted_species$taxon_name[which(wcvp_accepted_species$genus %in% unique_genera_in_tree)])
total_number_accepted_species

# Find the proportion of missing sp
missing_sp_biomes[1] <- 1-(missing_sp_biomes[1]/nr_trf)
missing_sp_biomes[2] <-1-(missing_sp_biomes[2]/nr_non_trf)
missing_sp_biomes[3] <- 1-(missing_sp_biomes[3]/nr_widespread)


cat("The proportion of missing species per biome are the following: \n")
cat("Trf: ", missing_sp_biomes[1], "\n")
cat("Non-trf: ", missing_sp_biomes[2], "\n")
cat("Widespread: ", missing_sp_biomes[3], "\n\n")

cat("Writing the output file \n\n")

# Now I need to save the result summary dataframe to a file which can be read by ESSE
write.table(missing_sp_biomes, file = paste0(path_out, output), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


################################################################################################################################################
######################################################-- End of script --#######################################################################
################################################################################################################################################


################################################################################################################################################
######################################################-- Outdated code --#######################################################################
################################################################################################################################################

  #missing_sp <- setdiff(wcvp_accepted_species_subset$taxon_name, tree$tip.label)
  
  #Finding the climate description for the missing species
  # missing_sp_climate <- wcvp_accepted_species_subset[which(wcvp_accepted_species_subset$taxon_name %in% missing_sp),]

  #Now we use our estimations of wide spread Trf and non-trf to estimate how these climate descriptions are distributed.
  # total_wet_tropical <- length(missing_sp_climate$climate_description[which(missing_sp_climate$climate_description == "wet tropical")])
  # wcvp_wet_tropical_widespread <- total_wet_tropical * Wet_tropical_widespread_percentage # The number of species with Wet tropical that are possible widespread
  # wcvp_wet_tropical_trf <- total_wet_tropical * Wet_tropical_trf_percentage # The number of species with Wet tropical that are possible Trf
  # wcvp_wet_tropical_non_trf <- total_wet_tropical * Wet_tropical_non_trf_percentage # The number of species with Wet tropical that are possible non-trf

  # total_non_wet_tropical <- length(missing_sp_climate$climate_description[which(missing_sp_climate$climate_description != "wet tropical")])
  # wcvp_non_wet_tropical_widespread <- total_non_wet_tropical * non_wet_tropical_widespread_percentage # The number of species with non-wet tropical that are possible widespread
  # wcvp_non_wet_tropical_trf <- total_non_wet_tropical * non_wet_tropical_trf_percentage # The number of species with non-wet tropical that are possible Trf
  # wcvp_non_wet_tropical_non_trf <- total_non_wet_tropical * non_wet_tropical_non_trf_percentage # The number of species with non-wet tropical that are possible non-trf

  # With these numbers we can now calculate the sampling fraction for Trf, non-trf and widespread
  
  # Find the number of tips that are Trf, non-trf and widespread
  # all(tree$tip.label %in% result_summary$wcvp_taxon_name) # TRUE

  # Find the accepted species which are not in the tree.
  #wcvp_accepted_species_subset <- wcvp_accepted_species[which(wcvp_accepted_species$taxon_name %in% tree$tip.label),]
  #species_not_in_tree <- setdiff(wcvp_accepted_species$taxon_name, tree$tip.label)
  

  # We need to find out how many of the species missing from the tree which are Trf, non-trf and widespread
  #tree_tips_not_in_occurrences <- tree$tip.label[which(!(tree$tip.label %in% result_summary$wcvp_taxon_name))]
  #all(tree_tips_not_in_occurrences %in% missing_sp_climate$taxon_name) # TRUE

  # We can now calculate probability of these tips with only wcvp data being Trf, non-trf and widespread
# for (species in tree_tips_not_in_occurrences) {
#   climate_description <- missing_sp_climate$climate_description[missing_sp_climate$taxon_name == species]
  
#   if (climate_description == "wet tropical") {
#     probabilities <- c(Wet_tropical_trf_percentage, Wet_tropical_non_trf_percentage, Wet_tropical_widespread_percentage)
#   } else {
#     probabilities <- c(non_wet_tropical_trf_percentage, non_wet_tropical_non_trf_percentage, non_wet_tropical_widespread_percentage)
#   }
  
#   selected_category <- sample(c("10", "01", "11"), size = 1, prob = probabilities)
  
#   if (selected_category == "10") {
#     result_summary <- rbind(result_summary, data.frame(wcvp_taxon_name = species, occurrences_trf = 1, occurrences_non_trf = 0, proportion_in_tropical_rainforest = 1, proportion_outside_tropical_rainforest = 0, climate_description = climate_description, trf = 1, non_trf = 0))
#   } else if (selected_category == "01") {
#     result_summary <- rbind(result_summary, data.frame(wcvp_taxon_name = species, occurrences_trf = 0, occurrences_non_trf = 1, proportion_in_tropical_rainforest = 0, proportion_outside_tropical_rainforest = 1, climate_description = climate_description, trf = 0, non_trf = 1))
#   } else {
#     result_summary <- rbind(result_summary, data.frame(wcvp_taxon_name = species, occurrences_trf = 1, occurrences_non_trf = 1, proportion_in_tropical_rainforest = 0.5, proportion_outside_tropical_rainforest = 0.5, climate_description = climate_description, trf = 1, non_trf = 1))
#   }
# }

# result_summary