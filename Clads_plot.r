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
	clads_tip_lambda <- rbind(clads_tip_lambda, data.frame(order = order_name, tip_label = tree$tip.label, lambda = CladsOutput$lambdatip_map))
	
	# print dim of the dataframe to keep track of the progress
	print(dim(clads_tip_lambda))

	# Remove the CladsOutput object from the environment
	rm(CladsOutput)
}

# Finding orders in the ClaDs output.
unique_orders <- unique(clads_tip_lambda$order)
unique_orders <- unique_orders[unique_orders != "Berberidopsidales"]
unique_orders <- unique_orders[unique_orders != "Buxales"]
unique_orders <- unique_orders[unique_orders != "Paracryphiales"]

unique_orders

# Specify the folder path
folder_path <- "/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data/"

# Get a list of all distribution files in the folder
distribution_files <- list.files(folder_path, pattern = ".*_distribution_data\\.txt", full.names = TRUE)


distribution_data <- data.frame()

# Loop through each file
for (i in 1:length(distribution_files)) {

	# Keeping track of the progress
	print(paste("Processing file", i, "of", length(distribution_files)))

	# Extract the file name
	file <- distribution_files[i]

	# Print order
	print(file)

	# Extracting the order from the file name
	order_name <- gsub("(.*)_distribution_data\\.txt", "\\1", basename(file))

	# Cheking if the order is in the ClaDs output
	if (order_name %in% unique_orders) {
		# Load the distribution file
		distribution_data_raw <- fread(file)
		
		# Append the distribution data to the distribution data frame
		distribution_data <- rbind(distribution_data, distribution_data_raw)
		
		# rm file
		rm(distribution_data_raw)
	} else {
		next
	}
}

# Changing "_" to " " in the tip_label column
clads_tip_lambda$tip_label <- gsub("_", " ", clads_tip_lambda$tip_label)

# Merging the 2 datasets
distribution_data_merged <- merge(clads_tip_lambda, distribution_data, by.x = "tip_label", by.y = "wcvp_taxon_name")
dim(distribution_data_merged)

head(clads_tip_lambda)
head(distribution_data)

# Sort the distribution_data_merged dataframe by lambda in descending order
distribution_data_merged_sorted <- distribution_data_merged[order(-distribution_data_merged$lambda), ]

# Print the first 100 values
head(distribution_data_merged_sorted, 100)

# Remove rows which have the order Metteniusales
distribution_data_merged_sorted <- distribution_data_merged_sorted[!distribution_data_merged_sorted$order == "Metteniusales", ]

distribution_data_merged <- distribution_data_merged_sorted

# Plot
# Add a column to the distribution data that indicates whether the proportion in tropical rainforest is greater than 70%
distribution_data_merged$in_tropical_rainforest60 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.6
distribution_data_merged$in_tropical_rainforest70 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.7
distribution_data_merged$in_tropical_rainforest80 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.8
distribution_data_merged$in_tropical_rainforest90 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.9

################################################################################################################################
############################----- Density functions of tip rate speciation in tropical rainforest -----#########################
################################################################################################################################

trfcol <- "chartreuse"
non_trfcol <- "deeppink"

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


# Density functions of tip rate speciation in tropical rainforest

densityplot60 <- ggplot(distribution_data_merged, aes(x = lambda, fill = in_tropical_rainforest60)) +
	geom_density(alpha = 0.5) +
	scale_fill_manual(values = c(non_trfcol, trfcol), labels = c("", "")) +
	ylab("Density") +
	xlab("Tip Rate Speciation") +
	theme_ipsum() +
	labs(title = "60%") +
	theme(legend.position = "none") +
	scale_x_log10()

densityplot70 <- ggplot(distribution_data_merged, aes(x = lambda, fill = in_tropical_rainforest70)) +
	geom_density(alpha = 0.5) +
	scale_fill_manual(values = c(non_trfcol, trfcol), labels = c("", "")) +
	ylab("Density") +
	xlab("Tip Rate Speciation") +
	theme_ipsum() +
	labs(title = "70%") +
	theme(legend.position = "none") +
	scale_x_log10()

densityplot80 <- ggplot(distribution_data_merged, aes(x = lambda, fill = in_tropical_rainforest80)) +
	geom_density(alpha = 0.5) +
	scale_fill_manual(values = c(non_trfcol, trfcol), labels = c("", "")) +
	ylab("Density") +
	xlab("Tip Rate Speciation") +
	theme_ipsum() +
	labs(title = "80%") +
	theme(legend.position = "none") +
	scale_x_log10()

densityplot90 <- ggplot(distribution_data_merged, aes(x = lambda, fill = in_tropical_rainforest90)) +
	geom_density(alpha = 0.5) +
	scale_fill_manual(values = c(non_trfcol, trfcol), labels = c("", "")) +
	ylab("Density") +
	xlab("Tip Rate Speciation") +
	theme_ipsum() +
	labs(title = "90%") +
	theme(legend.position = "none") +
	scale_x_log10()

# Combine all density plots into one grid
prow <- plot_grid(densityplot60, densityplot70, densityplot80, densityplot90, labels = c("A", "B", "C", "D"), label_size = 12, label_fontface = "bold")


#
# extract a legend that is laid out horizontally
legend_b <- get_legend(
  densityplot60 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

# add the legend underneath the row we made earlier. Give it 10%
# of the height of one plot (via rel_heights).
plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))
