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
##############################################-- Handling Command Line arguments --#############################################################
################################################################################################################################################

# Command line arguments
input_file_tree <- commandArgs(trailingOnly = TRUE)[1]
output <- commandArgs(trailingOnly = TRUE)[2]
input_file_wcvp <- commandArgs(trailingOnly = TRUE)[3]
path_out <- commandArgs(trailingOnly = TRUE)[4]
order <- commandArgs(trailingOnly = TRUE)[5]

################################################################################################################################################
######################################################-- Loading files --#######################################################################
################################################################################################################################################

# Load the wcvp dataset
wcvp <- readRDS(input_file_wcvp)
tree <- read.tree(input_file_tree)

# Filter the wcvp dataset to include only rows where taxon_status == "Accepted"
wcvp_accepted <- subset(wcvp, taxon_status == "Accepted")

################################################################################################################################################
###############################################-- Finding Environmental data --################################################################
################################################################################################################################################

# Create an empty dataframe to store the results
result_df <- data.frame(tip_name = character(), climate_description = character(), stringsAsFactors = FALSE)

# Loop through each tip in the tree
for (i in seq_along(tree$tip.label)) {
  # Get the species name from the tip
  species_name <- tree$tip.label[i]

  # Search for the species name in the wcvp_accepted dataset
  matching_row <- wcvp_accepted[wcvp_accepted$taxon_name == species_name, ]
  
  # If a match is found, record the climate description in the result dataframe
  if (nrow(matching_row) > 0) {
    climate_description <- matching_row$climate_description
    result_df <- rbind(result_df, data.frame(tip_name = tip, climate_description = climate_description, stringsAsFactors = FALSE))
  }
}

# Print the result dataframe
cat("The number of tips in the tree of ", order, " is ", length(tree$tip.label), "\n")


#if there are no NA's or Empty string in the climate column Convert Wet tropical to 1's and everything else but NA or "" to 0's
if(sum(is.na(result_df$climate_description)) == 0 & sum(result_df$climate_description == "") == 0){
  cat("There are no missing climate descriptions of the species in the tree")
  result_df$climate_description <- ifelse(result_df$climate_description == "Wet tropical", 1, 0)
  write.table(result_df, file.path(path_out, output), sep = "\t", row.names = FALSE)
} else {
	# Report how many species are lacking climate data
  cat("There are ", sum(is.na(result_df$climate_description)), " species with NA in the climate data \n")
  cat("There are ", sum(result_df$climate_description == ""), " species with empty string in the climate data \n")
  break
}
