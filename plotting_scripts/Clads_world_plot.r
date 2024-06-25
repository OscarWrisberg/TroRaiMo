
#Install necessary packages if not already installed
install.packages("ggplot2")
install.packages("dplyr")
install.packages("maps")
install.packages("gridExtra")

# Load libraries
library(ggplot2)
library(dplyr)
library(maps)
library(gridExtra)
library(data.table)

#################################################################################################################################################
# Example Data
set.seed(123)  # For reproducibility
latitudes <- rep(seq(-50, 70, by = 10), each = 40)  # Double the number of entries for each belt
lambda_values <- c(rnorm(260, mean = 1, sd = 0.5), rnorm(260, mean = 2, sd = 0.5))  # Adjusted number of values
biomes <- rep(c("Non-Trf", "Trf"), each = 260)  # Adjusted biomes 280
biomes <- sample(biomes)

# Create data frame
data <- data.frame(latitudes, lambda_values, biomes)
#################################################################################################################################################
# Now I need to gather my own data and prepare it for the plot.

# I need to use the code which loads all the tip rate speciation from all the trees.
# I then need to use the occurrences data in order to find out which latitudes each species is present in.

# Then I need to create a data frame with the following columns:
# 1. Latitude band
# 2. Lambda values
# 3. Biome
# 4. Species

# Based on this dataset I should then be able to use the below plotting script to create a nice violin plot with a world map background.
# which shows how the difference in tip rate speciation varies with latitude and biome.

#########################################################################################################################
############################################ Loading ClaDs runs on orders  ##############################################
#########################################################################################################################

# Loading the ClaDs run for the runs that ran on the Orders
# Specify the folder path
folder_path <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/orders" # Local
folder_path <- "/home/owrisberg/Trf_models/workflow/02_adding_orders/pruning/orders" # Srun

# Get a list of all Rdata files in the folder
file_list <- list.files(folder_path, pattern = "Clads_output_.*\\.Rdata", full.names = TRUE)
file_list

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
folder_path <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/subset_of_orders" # Local
folder_path <- "/home/owrisberg/Trf_models/workflow/02_adding_orders/pruning/subset_of_orders" # Srun

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

head(clads_tip_lambda)

###############################################################################################################################
############################################--- Biome affiliations ---#########################################################
###############################################################################################################################
# Finding orders in the ClaDs output.
unique_orders <- unique(clads_tip_lambda$order)
unique_orders <- unique_orders[unique_orders != "Berberidopsidales"]
unique_orders <- unique_orders[unique_orders != "Buxales"]
unique_orders <- unique_orders[unique_orders != "Paracryphiales"]


unique_orders

# Specify the folder path
folder_path <- "/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data/"
folder_path <- "/home/owrisberg/Trf_models/workflow/03_distribution_data" # Srun

# Get a list of all distribution files in the folder
distribution_files <- list.files(folder_path, pattern = ".*_distribution_data\\.txt", full.names = TRUE)

distribution_files

# Remove Typhaceae_distribution_data.txt
distribution_files <- distribution_files[!distribution_files == "/home/owrisberg/Trf_models/workflow/03_distribution_data/Typhaceae_distribution_data.txt"]

# Remove /home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data//Xyridaceae_Eriocaulaceae_distribution_data.txt"
distribution_files <- distribution_files[!distribution_files == "/home/owrisberg/Trf_models/workflow/03_distribution_data/Xyridaceae_Eriocaulaceae_distribution_data.txt"]

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

dim(distribution_data)
dim(clads_tip_lambda)

# Finding out which tips have distr values but no distribution data
clads_tip_lambda_no_distribution <- clads_tip_lambda[!clads_tip_lambda$tip_label %in% distribution_data$wcvp_taxon_name, ]
dim(clads_tip_lambda_no_distribution) # is there really 14 k tips with no distribution data?

# Find orders of the tips missing distribution data
unique(clads_tip_lambda_no_distribution$order

# Count how many species from each of these unqique orders are missing distribution data
missing_per_order <- data.frame(table(clads_tip_lambda_no_distribution$order))

# 


head(distribution_data)
head(clads_tip_lambda$tip_label)

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
distribution_data_merged$in_tropical_rainforest50 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.50
distribution_data_merged$in_tropical_rainforest60 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.66
distribution_data_merged$in_tropical_rainforest70 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.75
distribution_data_merged$in_tropical_rainforest80 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.80
distribution_data_merged$in_tropical_rainforest90 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.90

head(distribution_data_merged)

###############################################################################################################################
############################################--- Calculating latitudes for species ---#########################################################
###############################################################################################################################

# loading the renamed occurrence dataset.
occurrence_data <- readRDS("/home/au543206/GenomeDK/Trf_models/workflow/01_distribution_data/06_Renamed/gbif_renamed.rds") # this kills my pc

##################################################################################################################################################
# Get world map data
world_map <- map_data("world")

# Create the world map plot
map_plot <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "grey80", color = "grey80") +
  coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-60, 80)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "transparent", color = NA))

plot(map_plot)

# Create the violin plot
violin_plot <- ggplot(data, aes(x = lambda_values, y = factor(latitudes), fill = biomes)) +
	geom_violin(scale = "width", adjust = 1, alpha = 0.6) +
	scale_fill_manual(values = c("Trf" = "chartreuse", "Non-Trf" = "deeppink")) +
	theme_minimal() +
	theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank()) +  # Remove panel background
	labs(x = expression(lambda ~ "(lineages/myr)"), y = "Latitude") +
	theme(panel.background = element_rect(fill = "transparent", color = NA),
				axis.text.y = element_text(size = 10, color = "black"))

plot(violin_plot)

# Combine the plots
combined_plot <- ggplot() +
  # World map as background
  annotation_custom(ggplotGrob(map_plot), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  # Overlay violin plot
  annotation_custom(ggplotGrob(violin_plot), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

# Plot
print(combined_plot)
