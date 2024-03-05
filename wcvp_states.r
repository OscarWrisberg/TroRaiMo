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

# setwd("/home/owrisberg/Trf_models/workflow/02_adding_orders/pruning/orders") # srun
# input_file_tree <- "pruned_tree_order_Arecales_GBMB.tre"
# output <- "Test_arecales.txt"
# input_file_wcvp <- "/home/owrisberg/Trf_models/workflow/02_adding_orders/wcvp_names_apg_aligned.rds" #srun
# path_out <- "/home/owrisberg/Trf_models/workflow/03_distribution_data/" #srun
# order_in_question <- as.character("Arecales")
# apg  <- "../../../../TroRaiMo/apgweb_parsed.csv"
# renamed_occurence_file <- "/home/owrisberg/Trf_models/workflow/01_distribution_data/06_Renamed/gbif_renamed.rds" #srun
# koppen_biome_file <- "../../../../TroRaiMo/koppen_geiger_0p01.tif" #srun

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

# Print the command line arguments
cat("The input file for the tree is ", input_file_tree, "\n")
cat("The output file is ", output, "\n")
cat("The input file for the wcvp dataset is ", input_file_wcvp, "\n")
cat("The path to the output file is ", path_out, "\n")
cat("The order of the tree is ", order_in_question, "\n")


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

# # Find order function
# find_order <- function(fams, apg) {
#   fam_list <- character(0)
#   orders <- character(0)

#   for (family in unique(fams)) {
#     order <- apg[which(apg$Syn_Fam == family), "Clade"] # Finding the order of that family in APG file
#     order <- as.character(order[1]) # Selecting the order of the family
#     cat("family ", family, "Order", order, "\n") # Printing the family and order
#     fam_list <- c(fam_list, family)
#     orders <- c(orders, order)
#   }

#   df_orders <- data.frame(family = fam_list, order = orders)
#   return(df_orders)
# }

# # Running the function
# cat("Finding the order for each of the families \n")
# family_orders <- find_order(unique_families, apg)
# length(family_orders$order)

# # Merging the wcvp and family_orders data frames
# cat("Merging the wcvp and family_orders data frames \n")
# wcvp_accepted_species_orders <- merge(wcvp_accepted_species, family_orders, by.x = "family", by.y = "family")

# # Creating a subset of wcvp where we only have the accepted species in the order
# cat("Creating a subset of wcvp where we only have the accepted species in the order \n")
# wcvp_accepted_species_orders_subset <- wcvp_accepted_species_orders[which(wcvp_accepted_species_orders$order %in% order_in_question),]

#print(wcvp_accepted_species_orders_subset)

################################################################################################################################################
#############################################-- Finding Environmental data GBIF --##############################################################
################################################################################################################################################

# Subset the renamed_occurence dataset to only include the species in the tree
cat("Subsetting the renamed_occurence dataset to only include the species in the tree \n\n")
renamed_occurence_subset <- renamed_occurence[which(renamed_occurence$wcvp_taxon_name %in% tree$tip.label),]
dim(renamed_occurence_subset) 

#Remove occurences with NA in decimalLatitude or decimalLongitude
cat("Removing occurences with NA in decimalLatitude or decimalLongitude \n\n")
renamed_occurence_subset <- renamed_occurence_subset[which(!is.na(renamed_occurence_subset$decimalLatitude) & !is.na(renamed_occurence_subset$decimalLongitude)),]
dim(renamed_occurence_subset) # 76483

# How many of the species in the tree are found in the renamed_occurence_subset
cat("Out of ",length(tree$tip.label),"tips in the tree there are ",length(tree$tip.label[which(tree$tip.label %in% renamed_occurence_subset$wcvp_taxon_name)]),"of the species in the tree are found in the renamed_occurence_subset \n\n")

# Loading the Koppen biomes data
cat("Loading the Koppen biomes data \n\n")
koppen_biome_map <- rast(koppen_biome_file)

##################################################################################################################################################################
##################################################################################################################################################################

if (all(tree$tip.label %in% renamed_occurence_subset$wcvp_taxon_name)) { # If all the tip labels are found in the occurences use only them to create the presence absence matrix

  # Add a new column to the csv file which is Tropical rainforest or not.
  cat("Adding a new column to the csv file which is Tropical rainforest or not \n\n")
  koppen_biome_map_tropical_rainforest <- ifel(koppen_biome_map == 1 | koppen_biome_map == 2, 1, 0)

  # Extracting the biome for each of the occurrences
  renamed_occurence_subset$in_tropical_rainforest <- terra::extract(koppen_biome_map_tropical_rainforest, renamed_occurence_subset[,c("decimalLongitude", "decimalLatitude")])[,2]

  # Remove occurences with NA in in_tropical_rainforest
  cat("Removing ", length(which(is.na(terra::extract(koppen_biome_map_tropical_rainforest, renamed_occurence_subset[,c("decimalLongitude", "decimalLatitude")])[,2]))) ," occurences with NA in in_tropical_rainforest \n\n")
  renamed_occurence_subset <- renamed_occurence_subset[which(!is.na(renamed_occurence_subset$in_tropical_rainforest)),]

  # Summarize the number of occurrences which are inside and outside the tropical rainforest
  result_summary <- aggregate(renamed_occurence_subset$in_tropical_rainforest, by = list(renamed_occurence_subset$wcvp_taxon_name), FUN = function(x) c(sum(x == 1), sum(x == 0)))

  # Calculate the proportion of occurrences inside the tropical rainforest biome
  result_summary$proportion_in_tropical_rainforest <- result_summary$x[,1]/(result_summary$x[,1] + result_summary$x[,2])
  result_summary$proportion_outside_tropical_rainforest <- result_summary$x[,2]/(result_summary$x[,1] + result_summary$x[,2])

##################################################################################################################################################################
##################################################################################################################################################################

} else { # If all the species are not found in the occurrences, see if the species missing from the occurrences are found in the WCVP
  renamed_occurence_subset <- renamed_occurence[which(renamed_occurence$wcvp_taxon_name %in% tree$tip.label),]
  cat("Not all the species in the tree are found in the renamed_occurence_subset \n\n")
  cat("The following species are not found in the renamed_occurence_subset: \n")
  cat(setdiff(tree$tip.label, renamed_occurence_subset$wcvp_taxon_name), "\n\n")
  missing_sp <- setdiff(tree$tip.label, renamed_occurence_subset$wcvp_taxon_name)

  cat("Can we find a biome for the missing species in the wcvp dataset? \n\n")
  # subsetting the wcvp dataset to only include the species in the tree
  wcvp_subset <- wcvp[which(wcvp$taxon_name %in% tree$tip.label),]

  # Subsetting wcvp to only include accepted species and species with a climate description
  wcvp_subset <- wcvp_subset[which(wcvp_subset$taxon_status == "Accepted" & wcvp_subset$climate_description != "" & wcvp_subset$taxon_rank == "Species"),]
  cat("Are all the missing sp found in WCVP and do they have a Climate description: ", all(missing_sp %in% wcvp_subset$taxon_name), "\n")

##################################################################################################################################################################

  if (all(missing_sp %in% wcvp_subset$taxon_name)) { # If all the missing species are found in the wcvp use the wcvp climate column to add the missing species to the presence absence matrix

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
    result_summary <- aggregate(renamed_occurence_subset$in_tropical_rainforest, by = list(renamed_occurence_subset$wcvp_taxon_name), FUN = function(x) c(sum(x == 1), sum(x == 0)))

    # Splitting the second column, which is a column containing a dataframe with 2 columns, into two individual columns
    result_summary <- cbind(result_summary[,1] , as.data.frame(result_summary[,2]))
    
    # Changing the column names
    result_summary <- setNames(result_summary, c("wcvp_taxon_name", "occurrences_trf", "occurrences_non_trf"))

    # Calculate the proportion of occurrences inside the tropical rainforest biome
    result_summary$proportion_in_tropical_rainforest <- result_summary$occurrences_trf/(result_summary$occurrences_trf + result_summary$occurrences_non_trf)
    result_summary$proportion_outside_tropical_rainforest <- result_summary$occurrences_non_trf/(result_summary$occurrences_trf + result_summary$occurrences_non_trf)

    # Use the wcvp_subset to find the biome for the rest of the species
    cat("Using the wcvp_subset to find the biome for the rest of the species \n\n")

    # First I will find the species which are missing from the GBif dataset
    missing_sp <- setdiff(tree$tip.label, renamed_occurence_subset$wcvp_taxon_name)

    # I will then find their climate description in the wcvp_subset
    missing_sp_climate <- wcvp_subset[which(wcvp_subset$taxon_name %in% missing_sp),]

    # Now we convert the climate description from wcvp to an estimation of biome.
    missing_sp_climate$koppen_biome <- ifelse(missing_sp_climate$climate_description == "wet tropical", 1, 0)

    # If missing_sp_climate$koppen_biome is == 1 then the proportion_in_tropical_rainforest is 1, otherwise it is 0
    missing_sp_result_add <- data.frame(wcvp_taxon_name = missing_sp_climate$taxon_name, occurrences_trf = NA , occurrences_non_trf = NA, proportion_in_tropical_rainforest = ifelse(missing_sp_climate$koppen_biome == 1, 1, 0), proportion_outside_tropical_rainforest = ifelse(missing_sp_climate$koppen_biome == 1, 0, 1))


    # Now we can add the missing_sp_climate species to the result_summary
    result_summary <- rbind(result_summary, missing_sp_result_add)

##################################################################################################################################################################

  } else {
    cat("Not all the missing species are found in the wcvp dataset \n")
    cat("The following species are found in the wcvp dataset \n")
    cat(missing_sp, "\n")
    cat("The following species are not found in the wcvp dataset \n")
    cat(setdiff(missing_sp, wcvp_subset$taxon_name), "\n")
    wcvp_missing_sp <- setdiff(missing_sp, wcvp_subset$taxon_name)

    missing_sp_file_name <- paste0("Missing_sp_",order_in_question,".txt")

    cat("In order for the pipeline to progess were removing the missing species from the tree, but they can be found in a file called ", missing_sp_file_name, "\n")
    write.table(wcvp_missing_sp, file = paste0(path_out, "Missing_sp_",order_in_question,".txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

    # dropping the missing species from the tree.
    tree <- ape::drop.tip(tree, wcvp_missing_sp)

    # Starting over but this time with the missing species removed from the tree
    renamed_occurence_subset <- renamed_occurence[which(renamed_occurence$wcvp_taxon_name %in% tree$tip.label),]
    missing_sp <- setdiff(tree$tip.label, renamed_occurence_subset$wcvp_taxon_name)

    # subsetting the wcvp dataset to only include the species in the tree
    wcvp_subset <- wcvp[which(wcvp$taxon_name %in% tree$tip.label),]

    # Subsetting wcvp to only include accepted species and species with a climate description
    wcvp_subset <- wcvp_subset[which(wcvp_subset$taxon_status == "Accepted" & wcvp_subset$climate_description != "" & wcvp_subset$taxon_rank == "Species"),]
    cat("Are all the missing sp found in WCVP and do they have a Climate description: ", all(missing_sp %in% wcvp_subset$taxon_name), "\n")

####################################################################################################################################################################

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
    result_summary <- aggregate(renamed_occurence_subset$in_tropical_rainforest, by = list(renamed_occurence_subset$wcvp_taxon_name), FUN = function(x) c(sum(x == 1), sum(x == 0)))

    # Splitting the second column, which is a column containing a dataframe with 2 columns, into two individual columns
    result_summary <- cbind(result_summary[,1] , as.data.frame(result_summary[,2]))

    # Changing the column names
    result_summary <- setNames(result_summary, c("wcvp_taxon_name", "occurrences_trf", "occurrences_non_trf"))

    # Calculate the proportion of occurrences inside the tropical rainforest biome
    result_summary$proportion_in_tropical_rainforest <- result_summary$occurrences_trf/(result_summary$occurrences_trf + result_summary$occurrences_non_trf)
    result_summary$proportion_outside_tropical_rainforest <- result_summary$occurrences_non_trf/(result_summary$occurrences_trf + result_summary$occurrences_non_trf)

    # Use the wcvp_subset to find the biome for the rest of the species
    cat("Using the wcvp_subset to find the biome for the rest of the species \n\n")

    # First I will find the species which are missing from the GBif dataset
    missing_sp <- setdiff(tree$tip.label, renamed_occurence_subset$wcvp_taxon_name)

    # I will then find their climate description in the wcvp_subset
    missing_sp_climate <- wcvp_subset[which(wcvp_subset$taxon_name %in% missing_sp),]

    # Now we convert the climate description from wcvp to an estimation of biome.
    missing_sp_climate$koppen_biome <- ifelse(missing_sp_climate$climate_description == "wet tropical", 1, 0)

    # If missing_sp_climate$koppen_biome is == 1 then the proportion_in_tropical_rainforest is 1, otherwise it is 0
    missing_sp_result_add <- data.frame(wcvp_taxon_name = missing_sp_climate$taxon_name, occurrences_trf = NA , occurrences_non_trf = NA, proportion_in_tropical_rainforest = ifelse(missing_sp_climate$koppen_biome == 1, 1, 0), proportion_outside_tropical_rainforest = ifelse(missing_sp_climate$koppen_biome == 1, 0, 1))
    
    # Now we can add the missing_sp_climate species to the result_summary
    result_summary <- rbind(result_summary, missing_sp_result_add)

  }
}



# Now I need to save the result summary dataframe to a file which can be read by ESSE
write.table(result_summary, file = paste0(path_out, output), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

