################################################################################################################################################
######################################################-- Loading packages --####################################################################
################################################################################################################################################
# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("data.table", "ape", "phytools", "geiger", "castor", "MonoPhy", "terra", "dplyr", "ggplot2")

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

setwd("/home/owrisberg/Trf_models/workflow/02_adding_orders/pruning/orders") # srun
input_file_tree <- "pruned_tree_order_Arecales_GBMB.tre"
output <- "Test_arecales.txt"
input_file_wcvp <- "/home/owrisberg/Trf_models/workflow/02_adding_orders/wcvp_names_apg_aligned.rds" #srun
path_out <- "/home/owrisberg/Trf_models/workflow/03_distribution_data/" #srun
order_in_question <- as.character("Arecales")
apg  <- "../../../../TroRaiMo/apgweb_parsed.csv"
renamed_occurence_file <- "/home/owrisberg/Trf_models/workflow/01_distribution_data/05_Taxon_match/gbif_renamed.rds" #srun
koppen_biome_file <- "../../../../TroRaiMo/koppen_geiger_0p01.tif" #srun

################################################################################################################################################
##############################################-- Handling Command Line arguments --#############################################################
################################################################################################################################################

# Command line arguments
# input_file_tree <- commandArgs(trailingOnly = TRUE)[1]
# output <- commandArgs(trailingOnly = TRUE)[2]
# input_file_wcvp <- commandArgs(trailingOnly = TRUE)[3]
# path_out <- commandArgs(trailingOnly = TRUE)[4]
# order_in_question <- commandArgs(trailingOnly = TRUE)[5]
# apg <- commandArgs(trailingOnly = TRUE)[6]
# renamed_occurence <- commandArgs(trailingOnly = TRUE)[7]
# koppen_biome <- commandArgs(trailingOnly = TRUE)[8]

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

# Find order function
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
cat("Finding the order for each of the families \n")
family_orders <- find_order(unique_families, apg)
length(family_orders$order)

# Merging the wcvp and family_orders data frames
cat("Merging the wcvp and family_orders data frames \n")
wcvp_accepted_species_orders <- merge(wcvp_accepted_species, family_orders, by.x = "family", by.y = "family")

# Creating a subset of wcvp where we only have the accepted species in the order
cat("Creating a subset of wcvp where we only have the accepted species in the order \n")
wcvp_accepted_species_orders_subset <- wcvp_accepted_species_orders[which(wcvp_accepted_species_orders$order %in% order_in_question),]

print(wcvp_accepted_species_orders_subset)

################################################################################################################################################
#############################################-- Finding Environmental data GBIF --##############################################################
################################################################################################################################################

# Subset the renamed_occurence dataset to only include the species in the tree
cat("Subsetting the renamed_occurence dataset to only include the species in the tree \n\n")
renamed_occurence_subset <- renamed_occurence[which(renamed_occurence$wcvp_taxon_name %in% tree$tip.label),]


# Loading the Koppen biomes data
cat("Loading the Koppen biomes data \n\n")
koppen_biome_map <- rast(koppen_biome_file)

# Add a new column to the csv file which is Tropical rainforest or not.
cat("Adding a new column to the csv file which is Tropical rainforest or not \n\n")
koppen_biome_map_tropical_rainforest <- ifel(koppen_biome_map == 1 | koppen_biome_map == 2, 1, 0)


# Converting the latitudes and longitudes of the occurrences to spatial points
cat("Converting the latitudes and longitudes to spatial points \n\n")
coordinates <- c("decimalLongitude", "decimalLatitude")

# Create a SpatVector object and defining crs
spatial_points <- vect(renamed_occurence_subset, geom = coordinates)
crs(spatial_points) <- "+proj=longlat +datum=WGS84"

# Extracting the koppen biomes for each of the occurences
cat("Extracting the koppen biomes for each of the occurences \n\n")
extracted_biomes <- extract(koppen_biome_map_tropical_rainforest, spatial_points)
renamed_occurence_subset$koppen_biome <- extracted_biomes

# Counting the number of occurrences in each biome per species
cat("Counting the number of occurrences in each biome per species \n\n")
biome_count <- renamed_occurence_subset[, .(count = .N), by = .(wcvp_taxon_name, koppen_biome)]

points_df <- as.data.frame(spatial_points)

# Create a ggplot using the raster data
plot1 <- ggplot() +
  geom_raster(data = terra::vect(koppen_biome_map_tropical_rainforest), aes(x = x, y = y, fill = layer)) +
  
  # Add spatial points
  geom_point(data = points_df, aes(x = x, y = y), color = "red", size = 3) +
  
  # Customize plot appearance
  scale_fill_gradient(low = "white", high = "blue") +  # Adjust color scale for raster
  theme_minimal()

pdf("koppen_biome_map.pdf")
plot(koppen_biome_map_tropical_rainforest)
plot(spatial_points, col="red", pch=16, cex=0.5)
dev.off()

# pdf("spatial_points.pdf")
# plot(spatial_points, col="red", pch=16, cex=1)
# dev.off()


if (all(tree$tip.label %in% renamed_occurence_subset$wcvp_taxon_name)) {
  
# THings that I probably need to do in this step:
# Convert the latitudes and longitudes to spatial points.
# check for "bad" spatial points IE. points that are not on land. ( could potentially use a R package for this)
# Calculate the number of occurrences falling inside biomes 12 and 13 (tropical rainforest) and the number of occurrences falling outside biomes 12 and 13 (not tropical rainforest).
# Calculate the proportion of occurrences falling inside biomes 12 and 13 (tropical rainforest) and the proportion of occurrences falling outside biomes 12 and 13 (not tropical rainforest).
# This could then be used to create some different files with the in which we have the presence-absence of the species in the different biomes.
  # We would make this presence absence file for several different cut-off points for the proportion of species occurrences needing to be in a biome in order for it to be registered as present.
  # These cutoffs could be 0.1, 0.25

  # Loading the Koppen biomes data
  cat("Loading the Koppen biomes data \n\n")
  koppen_biome_map <- read.csv(koppen_biome_file, sep = ",", header = TRUE)



} else {
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
  head(wcvp_subset)
  cat("Are all the missing sp found in WCVP and do they have a Climate description: ", all(missing_sp %in% wcvp_subset$taxon_name), "\n")
  if (all(missing_sp %in% wcvp_subset$taxon_name)) {
    # use the koppen biome maps for the species represented in the renamed_occurence_subset

    # Use the wcvp_subset to find the biome for the rest of the species
  } else {
    cat("Not all the missing species are found in the wcvp dataset \n")
    cat("The following species are found in the wcvp dataset \n")
    cat(missing_sp, "\n")
    cat("The following species are not found in the wcvp dataset \n")
    cat(setdiff(missing_sp, wcvp_subset$taxon_name), "\n")
    break
  }
}





species_missing_in_gbif <- tree$tip.label[which(!tree$tip.label %in% renamed_occurence_subset$wcvp_taxon_name)]


# Loading the Koppen biomes data
cat("Loading the Koppen biomes data \n\n")
koppen_biome_map <- read.csv(koppen_biome, sep = ",", header = TRUE)

# Add a new column to the csv file which is Tropical rainforest or not.
cat("Adding a new column to the csv file which is Tropical rainforest or not \n\n")
koppen_biome_map$koppen <- ifelse(koppen_biome_map$koppen == 12 | koppen_biome_map$koppen == 13, 1, 0)

# Converting the latitudes and longitudes to spatial points
cat("Converting the latitudes and longitudes to spatial points \n\n")
coordinates(renamed_occurence_subset)  <- c(decimalLongitude,decimalLatitude) 


################################################################################################################################################
#############################################-- Finding Environmental data WCVP --##############################################################
################################################################################################################################################

# Create an empty dataframe to store the results
result_df <- data.frame(taxon_name = character(), climate_description = character(), stringsAsFactors = FALSE)

# Check if all the tips in the tree are found in the wcvp_accepted_species_orders
cat("Are all the tips in the tree found in the wcvp_accepted_species_orders ", all(tree$tip.label %in% wcvp_accepted_species_orders$taxon_name), "\n")
if(all(tree$tip.label %in% wcvp_accepted_species_orders$taxon_name) == FALSE){
  cat("Not all the tips in the tree are found in the wcvp_accepted_species_orders \n")
  cat("The following tips are not found in the wcvp_accepted_species_orders \n")
  cat(setdiff(tree$tip.label, wcvp_accepted_species_orders$taxon_name), "\n")
}


# Loop through each tip in the wcvp_accepted_species_orders
for (i in seq_along(wcvp_accepted_species_orders_subset$taxon_name)) {
  # Get the species name from the tip
  species_name <- wcvp_accepted_species_orders_subset$taxon_name[i]
  #cat("The species name is ", species_name, "\n")

  # Search for the species name in the wcvp_accepted dataset
  matching_row <- wcvp_accepted[wcvp_accepted$taxon_name == species_name, ]
  #cat("This is the matching row: \n")
  #print(matching_row)

  # If a match is found, record the climate description in the result dataframe
  if (nrow(matching_row) > 0) {
    climate_description <- matching_row$climate_description
    result_df <- rbind(result_df, data.frame(taxon_name = species_name, climate_description = climate_description, stringsAsFactors = FALSE))
  }
}

# Print the result dataframe
cat("The number of tips in the tree of ", order_in_question, " is ", length(tree$tip.label), "\n")


#if there are no NA's or Empty string in the climate column Convert Wet tropical to 1's and everything else but NA or "" to 0's
if(sum(is.na(result_df$climate_description)) == 0 & sum(result_df$climate_description == "") == 0){
  cat("There are no missing climate descriptions of the species in the wcvp dataset \n")
  result_df$climate_description <- ifelse(result_df$climate_description == "wet tropical", 1, 0)
  write.table(result_df, file.path(path_out, output), sep = "\t", row.names = FALSE)
} else {
	# Report how many species are lacking climate data
  cat("There are ", sum(is.na(result_df$climate_description)), " species with NA in the climate data \n")
  cat("There are ", sum(result_df$climate_description == ""), " species with empty string in the climate data \n")
  break
}


