# In this file I want to loop through all the ClaDs output files and see if the reason for finishing is because they ran out of time.

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

########################################################################################################################
############################################## Checking for duplicates  ################################################
########################################################################################################################

# Clads on orders
folderpath_1 <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/orders"
file_list_1 <- list.files(folderpath_1, pattern = "Clads_output_.*\\.Rdata", full.names = TRUE)

#Clads on suborders
folderpath_2 <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/subset_of_orders"
file_list_2 <- list.files(folderpath_2, pattern = "Clads_output_.*\\.Rdata", full.names = TRUE)

# join the file lists
file_list <- c(file_list_1, file_list_2)

file_list # length = 264 30/08/2024

########################################################################################################################
############################################## Creating the dataframe  ################################################
########################################################################################################################
# Loop through each file in the folder
for (i in 1:length(file_list)) {

	# Keeping track of the progress
	print(paste("Processing file", i, "of", length(file_list)))

	# Extractinc the file name
	file = file_list[i]

	# Extract the order name from the file name
	order_name <- gsub("Clads_output_(.*)\\.Rdata", "\\1", basename(file))
	
	# Load the Rdata file
	load(file)

	if ("reason_for_stop" %in% names(CladsOutput)) {
		#print("The reason for stop is in the file")
		if(CladsOutput$reason_for_stop == "max_time") {
			print("The reason for stop is because the run ran out of time")
			Unfinished_runs <- c(Unfinished_runs, order_name)
		}
	} else {
		#print("The reason for stop is not in the file")
		next
	}

	# Remove the CladsOutput object from the environment
	rm(CladsOutput)
}
