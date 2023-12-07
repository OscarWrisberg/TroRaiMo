################################################################################################################################################
######################################################-- Loading packages --####################################################################
################################################################################################################################################
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
################################################################################################################################################
#######################################################-- Local testing --######################################################################
################################################################################################################################################

# setwd("/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/")
# input_file_tree <- "pruned_tree__order_Arecales_GBMB.txt"
# output <- "Test_arecales.txt"
# input_file_wcvp <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/wcvp_names_apg_aligned.rds"
# path_out <- "/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data/"
# order_in_question <- as.character("Arecales")
# apg <- apg <- fread("../../../TroRaiMo/apgweb_parsed.csv")


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
###############################################-- Finding Environmental data --################################################################
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


