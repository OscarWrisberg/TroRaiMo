#########################################################################################################################
################################################## Loading packages  ####################################################
#########################################################################################################################

# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("data.table", "ape", "phytools", "geiger", "castor", "dplyr", "ggplot2", "viridis","hrbrthemes", "cowplot", "MetBrewer")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

#########################################################################################################################
################################################# Loading ClaDs runs  ###################################################
#########################################################################################################################


# Specify the folder path
folder_path <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/orders"

# Get a list of all Rdata files in the folder
file_list <- list.files(folder_path, pattern = "Clads_output_.*\\.Rdata", full.names = TRUE)

# Create dataframe for tips and tip rate speciation
clads_tip_lambda <- data.frame()

  
# Loop through each file in the folder
for (i in 1:length(file_list)) {

	# Keeping track of the progress
	print(paste("Processing file", i, "of", length(file_list)))

	# Extractinc the file name
	file = file_list[1]

	# Extract the order name from the file name
	order_name <- gsub("Clads_output_(.*)\\.Rdata", "\\1", basename(file))
	
	# Load the Rdata file
	load(file)

	tree <- CladsOutput$tree

	# Append the tip names and tip rate speciation to the dataframe
	clads_tip_lambda <- rbind(clads_tip_lambda, data.frame(order = order_name, tip_label = tree$tip.label, lambda = CladsOutput$lambdatip_map))
	
	# Remove the CladsOutput object from the environment
	rm(CladsOutput)
}


# Specify the folder path
folder_path <- "/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data/"

# Get a list of all distribution files in the folder
distribution_files <- list.files(folder_path, pattern = ".*_distribution_data\\.txt", full.names = TRUE)

distribution_data <- data.frame()

# Loop through each file
for (file in distribution_files) {
	print(file)
	# Load the distribution file
	distribution_data_raw <- fread(file)
	
	# Append the distribution data to the distribution data frame
	distribution_data <- rbind(distribution_data, distribution_data_raw)
	
}
distribution_data



distribution_data

# Filter the distribution data based on the proportion_in_tropical_rainforest column

# Create the density plot
ggplot(distribution_data, aes(x = proportion_in_tropical_rainforest)) +
	geom_density() +
	xlab("Proportion in Tropical Rainforest") +
	ylab("Density")


# Plot
# Add a column to the distribution data that indicates whether the proportion in tropical rainforest is greater than 70%
distribution_data_merged$in_tropical_rainforest60 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.6
distribution_data_merged$in_tropical_rainforest70 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.7
distribution_data_merged$in_tropical_rainforest80 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.8
distribution_data_merged$in_tropical_rainforest90 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.9

################################################################################################################################
############################----- Density functions of tip rate speciation in tropical rainforest -----#########################
################################################################################################################################

trfcol <- met.brewer("Cassatt2",7)[7]
non_trfcol <- met.brewer("Cassatt2",7)[2]

boxplot60 <- ggplot(distribution_data_merged, aes(x = in_tropical_rainforest60, y = lambda)) +
	geom_boxplot(fill = c(non_trfcol, trfcol), color = "black") +
	scale_x_discrete(labels = c("Not in TRF", "In TRF")) +
	ylab("Tip Rate Speciation") +
	xlab("") +
	theme_ipsum() + 
	labs(title = "60%") +
	theme(legend.position = "none")

boxplot70 <- ggplot(distribution_data_merged, aes(x = in_tropical_rainforest70, y = lambda)) +
	geom_boxplot(fill = c(non_trfcol, trfcol), color = "black") +
	scale_x_discrete(labels = c("Not in TRF", "In TRF")) +
	ylab("Tip Rate Speciation") +
	xlab("") +
	theme_ipsum() +
	labs(title = "70%") +
	theme(legend.position = "none")

boxplot80 <- ggplot(distribution_data_merged, aes(x = in_tropical_rainforest80, y = lambda)) +
	geom_boxplot(fill = c(non_trfcol, trfcol), color = "black") +
	scale_x_discrete(labels = c("Not in TRF", "In TRF")) +
	ylab("Tip Rate Speciation") +
	xlab("") +
	theme_ipsum() +
	labs(title = "80%") +
	theme(legend.position = "none")

boxplot90 <- ggplot(distribution_data_merged, aes(x = in_tropical_rainforest90, y = lambda)) +
	geom_boxplot(fill = c(non_trfcol, trfcol), color = "black") +
	scale_x_discrete(labels = c("Not in TRF", "In TRF")) +
	ylab("Tip Rate Speciation") +
	xlab("") +
	theme_ipsum() +
	labs(title = "90%") +
	theme(legend.position = "none")

# Combine all boxplots into one grid
plot_grid(boxplot60, boxplot70, boxplot80, boxplot90, labels = c("A", "B", "C", "D"), label_size = 12, label_fontface = "bold")

