#########################################################################################################################
################################################## Loading packages  ####################################################
#########################################################################################################################

# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("data.table", "ape", "phytools", "geiger", "castor", "dplyr")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

#########################################################################################################################
################################################# Loading ClaDs runs  ###################################################
#########################################################################################################################"

# Command line arguments
folder_path <- commandArgs(trailingOnly = TRUE)[1]
output_file <- commandArgs(trailingOnly = TRUE)[2]

# Get a list of all Rdata files in the folder
file_list <- list.files(folder_path, pattern = "Clads_output_.*\\.Rdata", full.names = TRUE)

# Create dataframe for tips and tip rate speciation
clads_tip_lambda <- data.frame()

# Loop through each file in the folder
for (i in 1:length(file_list)) { # This loop takes atleast 30 mins to run ....

	# Keeping track of the progress
	print(paste("Processing file", i, "of", length(file_list)))

	# Extracting the file name
	file <- file_list[i]

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
	
	# Print dim of the dataframe to keep track of the progress
	cat("The dataset contains ", dim(clads_tip_lambda)[1], " rows and ", dim(clads_tip_lambda)[2], " columns\n")

	# Remove the CladsOutput object from the environment
}

rm(CladsOutput)

# Extract the unique extinction rates
unique_extinction <- unique(clads_tip_lambda$extinction)

# Calculate the mean and standard deviation of the unique extinction rates
mean_extinction <- mean(unique_extinction)
log(mean_extinction)
sd_extinction <- sd(unique_extinction)
log(sd_extinction)

# Save the mean and sd extinction to a dataframe
extinction_stats <- data.frame(mean = mean_extinction, sd = sd_extinction)

# Save the dataframe to a csv file
write.csv(extinction_stats, file = output_file, row.names = FALSE)

