# Test 
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("maps")
# install.packages("gridExtra")

# Load libraries
library(ggplot2)
library(dplyr)
library(maps)
library(gridExtra)
library(data.table)
library(cowplot)
library(dplyr)

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

saveRDS(clads_tip_lambda, "/home/owrisberg/Trf_models/workflow/01_distribution_data/06_Renamed/clads_tip_lambda.rds") # Saving the clads_tip_lambda dataframe when working on the server
clads_tip_lambda <- readRDS("/home/au543206/GenomeDK/Trf_models/workflow/01_distribution_data/06_Renamed/clads_tip_lambda.rds") # loading it for working locally.
clads_tip_lambda <- readRDS("/home/owrisberg/Trf_models/workflow/01_distribution_data/06_Renamed/clads_tip_lambda.rds")

# Remove the _ from the tip_label column
clads_tip_lambda$tip_label <- gsub("_", " ", clads_tip_lambda$tip_label)

###############################################################################################################################
############################################--- Biome affiliations ---#########################################################
###############################################################################################################################
# Finding orders in the ClaDs output.
unique_orders <- unique(clads_tip_lambda$order)
unique_orders <- unique_orders[unique_orders != "Berberidopsidales"]
unique_orders <- unique_orders[unique_orders != "Buxales"]
unique_orders <- unique_orders[unique_orders != "Paracryphiales"]


unique_orders

# in unique_orders remove sub_phylo and .tre
unique_orders <- gsub("sub_phylo_", "", unique_orders)
unique_orders <- gsub(".tre", "", unique_orders)

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

# Remove Asteraceae_13_distribution_data.txt
distribution_files <- distribution_files[!distribution_files == "/home/owrisberg/Trf_models/workflow/03_distribution_data/Asteraceae_13_distribution_data.txt"]

distribution_data <- data.frame()
column_names <- c("wcvp_taxon_name", "occurrences_trf", "occurrences_non_trf", "proportion_in_tropical_rainforest", "proportion_outside_tropical_rainforest")

names(distribution_data) <- column_names

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

		# if the column names are not 
		if (!all(names(distribution_data_raw) %in% column_names)) {
			names(distribution_data_raw) <- column_names
		}
		
		# Append the distribution data to the distribution data frame
		distribution_data <- rbind(distribution_data, distribution_data_raw)
		
		# rm file
		rm(distribution_data_raw)
	} else {
		next
	}
}


dim(distribution_data)
dim(clads_tip_lambda)

# Finding out which tips have distr values but no distribution data
clads_tip_lambda_no_distribution <- clads_tip_lambda[!clads_tip_lambda$tip_label %in% distribution_data$wcvp_taxon_name, ]
dim(clads_tip_lambda_no_distribution) # 1206 tips with no distribution

# Find orders of the tips missing distribution data
unique(clads_tip_lambda_no_distribution$order)

# Count how many species from each of these unqique orders are missing distribution data
missing_per_order <- data.frame(table(clads_tip_lambda_no_distribution$order))
missing_per_order

head(distribution_data)
head(clads_tip_lambda$tip_label)

# Merging the 2 datasets
distribution_data_merged <- merge(clads_tip_lambda, distribution_data, by.x = "tip_label", by.y = "wcvp_taxon_name")
dim(distribution_data_merged)

head(clads_tip_lambda)
head(distribution_data)

dim(distribution_data_merged)

# Sort the distribution_data_merged dataframe by lambda in descending order
distribution_data_merged_sorted <- distribution_data_merged[order(-distribution_data_merged$lambda), ]

# Print the first 100 values
head(distribution_data_merged_sorted, 100)

# Remove rows which have the order Metteniusales
distribution_data_merged_sorted <- distribution_data_merged_sorted[!distribution_data_merged_sorted$order == "Metteniusales", ]

distribution_data_merged <- distribution_data_merged_sorted

# save this distribution_data_merged for later use
saveRDS(distribution_data_merged, "/home/owrisberg/Trf_models/workflow/01_distribution_data/06_Renamed/distribution_data_merged.rds")
distribution_data_merged <- readRDS("/home/au543206/GenomeDK/Trf_models/workflow/01_distribution_data/06_Renamed/distribution_data_merged.rds")
distribution_data_merged <- readRDS("/home/owrisberg/Trf_models/workflow/01_distribution_data/06_Renamed/distribution_data_merged.rds")

# Add a column to the distribution data that indicates whether the proportion in tropical rainforest is greater than 70%
distribution_data_merged$in_tropical_rainforest50 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.50
distribution_data_merged$in_tropical_rainforest60 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.66
distribution_data_merged$in_tropical_rainforest70 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.75
distribution_data_merged$in_tropical_rainforest80 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.80
distribution_data_merged$in_tropical_rainforest90 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.90

# Add a column showing whether the species is in tropical rainforest at 0.33
distribution_data_merged$in_tropical_rainforest33 <- distribution_data_merged$proportion_in_tropical_rainforest > 0.33

# Add a column showing wether the species is outside tropical rainforest at 0.33
distribution_data_merged$outside_tropical_rainforest33 <- distribution_data_merged$proportion_outside_tropical_rainforest > 0.33

# Now we add a column which shows either Trf, Non-trf or widespread based on the values of in_tropical_rainforest33 and outside_tropical_rainforest33
distribution_data_merged$biome_affiliation <- ifelse(
  distribution_data_merged$in_tropical_rainforest33 == TRUE & distribution_data_merged$outside_tropical_rainforest33 == FALSE,
  "TRF",
  ifelse(
    distribution_data_merged$outside_tropical_rainforest33 == TRUE & distribution_data_merged$in_tropical_rainforest33 == FALSE,
    "Non TRF",
    ifelse(
      distribution_data_merged$in_tropical_rainforest33 == TRUE & distribution_data_merged$outside_tropical_rainforest33 == TRUE,
      "Widespread",
      "Error"
    )
  )
)

###############################################################################################################################
############################################--- Calculating latitudes for species ---#########################################################
###############################################################################################################################

# Loading the renamed occurrence dataset.
occurrence_data <- readRDS("/home/owrisberg/Trf_models/workflow/01_distribution_data/06_Renamed/gbif_renamed.rds") # this kills my pc


# Changing "_" to " " in the tip_label column
clads_tip_lambda$tip_label <- gsub("_", " ", clads_tip_lambda$tip_label)

# create subset which only includes tips that are in the clads_tip_lambda dataframe
occurrence_data_subset <- occurrence_data[occurrence_data$wcvp_taxon_name %in% clads_tip_lambda$tip_label, ]
dim(occurrence_data_subset)

# Now we loop through the latitudinal bands and find out which species are present in each band
# We then make dataframe for each of these latitudinal bands which includes the species and their lambda values, the latitude band and the biome affiliation.

# Create a dataframe to store the data
data_for_plot <- data.frame()

# Loop through the latitudinal bands
latitudinal_bands <- seq(-50, 70, by = 10)
for (i in 1:(length(latitudinal_bands) - 1)) {
	
	# Keeping track of the progress
	print(paste("Processing latitude band", i, "of", length(latitudinal_bands) - 1))

	# Extract the lower and upper bounds of the latitude band
	lower_bound <- latitudinal_bands[i]
	print(lower_bound)
	upper_bound <- latitudinal_bands[i + 1]
	print(upper_bound)

	# Extract the species in the latitude band
	occurrence_data_subset_lat_band <- occurrence_data_subset[occurrence_data_subset$decimalLatitude >= lower_bound & occurrence_data_subset$decimalLatitude < upper_bound,]

	# Keep only unique species from the latitude band
	occurrence_data_subset_lat_band <- unique(occurrence_data_subset_lat_band$wcvp_taxon_name)

	# subset the Clads data based on the occurrence data
	occurrence_data_subset_lat_band_clads <- clads_tip_lambda[clads_tip_lambda$tip_label %in% occurrence_data_subset_lat_band,]

	head(occurrence_data_subset_lat_band_clads)

	# Append the data to the data_for_plot dataframe
	data_for_plot <- rbind(data_for_plot, data.frame(lat_band = paste(lower_bound, upper_bound, sep = "_"),
													  wcvp_taxon_name = occurrence_data_subset_lat_band_clads$tip_label,
													  lambda = occurrence_data_subset_lat_band_clads$lambda,
													  order = occurrence_data_subset_lat_band_clads$order))
}


saveRDS(data_for_plot, "/home/owrisberg/Trf_models/workflow/01_distribution_data/06_Renamed/data_for_plot.rds") # Saving the data_for_plot dataframe when working on the server
data_for_plot <- readRDS("/home/au543206/GenomeDK/Trf_models/workflow/01_distribution_data/06_Renamed/data_for_plot.rds") # loading it for working locally.
data_for_plot <- readRDS("/home/owrisberg/Trf_models/workflow/01_distribution_data/06_Renamed/data_for_plot.rds")

head(data_for_plot)
head(distribution_data_merged)

# Now I add the biome affiliation to the data_for_plot dataframe
data_for_plot <- merge(data_for_plot, distribution_data_merged, by.x = "wcvp_taxon_name", by.y = "tip_label")

head(data_for_plot)

##################################################################################################################################################

# Selecting colours for the plot
cols <- met.brewer("Java", n=3, type = "discrete")

#change the order of the colours
cols <- cols[c(1,3,2)]


data_for_plot_counts <- data.frame()

# Calculate counts for each group and for each level of trf affiliation
data_for_plot_counts_60 <- data_for_plot %>%
	group_by(lat_band, in_tropical_rainforest60) %>%
	summarize(count = n())

data_for_plot_counts_70 <- data_for_plot %>%
	group_by(lat_band, in_tropical_rainforest70) %>%
	summarize(count = n())

data_for_plot_counts_80 <- data_for_plot %>%
	group_by(lat_band, in_tropical_rainforest80) %>%
	summarize(count = n())

data_for_plot_counts_90 <- data_for_plot %>%
	group_by(lat_band, in_tropical_rainforest90) %>%
	summarize(count = n())

data_for_plot_counts_biome_affiliation <- data_for_plot %>%
	group_by(lat_band, biome_affiliation) %>%
	summarize(count = n())


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

unique(data_for_plot$lat_band)

# Define the desired order of lat_band levels
desired_order <- c("-50_-40","-40_-30","-30_-20","-20_-10","-10_0","0_10","10_20","20_30","30_40","40_50","50_60","60_70")  # Replace with your actual levels

# Convert lat_band to a factor with the specified levels
data_for_plot$lat_band <- factor(data_for_plot$lat_band, levels = desired_order)


# # Create the violin plot
# violin_plot <- ggplot(data_for_plot, aes(x = log10(lambda.x), y = lat_band, fill = in_tropical_rainforest90)) +
# 	geom_violin(scale = "width", adjust = 1, alpha = 0.6) +
# 	scale_fill_manual(values = c("TRUE" = "chartreuse", "FALSE" = "deeppink")) +
# 	theme_minimal() +
# 	theme(
#     panel.grid.major = element_blank(),  # Remove major grid lines
#     panel.grid.minor = element_blank(),  # Remove minor grid lines
#     panel.background = element_blank()) +  # Remove panel background
# 	labs(x = expression(lambda ~ "(lineages/myr)"), y = "Latitude") +
# 	theme(panel.background = element_rect(fill = "transparent", color = NA), axis.text.y = element_text(size = 10, color = "black"))

# plot(violin_plot)

# # Create the boxplot
# boxplot <- ggplot(data_for_plot, aes(x = log(lambda.x), y = lat_band, fill = as.factor(in_tropical_rainforest90))) +
#   geom_boxplot(alpha = 0.6) +
#   scale_fill_manual(values = c("TRUE" = "chartreuse", "FALSE" = "deeppink"),
#   labels = c("TRUE" = "TRF", "FALSE" = "Non TRF")) +
#   theme_minimal() +
#   theme(
#     panel.grid.major = element_blank(),  # Remove major grid lines
#     panel.grid.minor = element_blank(),  # Remove minor grid lines
#     panel.background = element_blank()  # Remove panel background
#   ) +
#   labs(x = expression(lambda ~ "log((lineages/myr)")), y = "Latitude") +
#   theme(panel.background = element_rect(fill = "transparent", color = NA), axis.text.y = element_text(size = 10, color = "black"))
    
# plot(boxplot)


# # Create the boxplot
# boxplot <- ggplot(data_for_plot, aes(x = log(lambda.x), y = lat_band, fill = as.factor(in_tropical_rainforest70))) +
#   geom_boxplot(alpha = 0.6) +
#   scale_fill_manual(values = c("TRUE" = "chartreuse", "FALSE" = "deeppink"),
#   					labels = c("TRUE" = "TRF", "FALSE" = "Non TRF")) +
#   theme_minimal() +
#   theme(
#     panel.grid.major = element_blank(),  # Remove major grid lines
#     panel.grid.minor = element_blank(),  # Remove minor grid lines
#     panel.background = element_blank()  # Remove panel background
#   ) +
#   labs(x = expression(lambda ~ "log((lineages/myr)"), y = "Latitude") +
#   #ggtitle("Tip rate speciation across latitudinal bands") +
#   theme(
#     legend.position = "bottom",
#     legend.title = element_blank(),
#     legend.text = element_text(size = 10),
#     legend.key = element_rect(fill = NA, color = NA),
#     legend.key.size = unit(1.5, "lines"),
#     panel.background = element_rect(fill = "transparent", color = NA),
#     axis.text.y = element_text(size = 10, color = "black")
#   ) +
#   geom_text(data = data_for_plot_counts_70, aes(label = count, x = 9, y = lat_band),
#             position = position_dodge(width = 0.75), vjust = -0.5, color = "black")

# plot(boxplot)


# Find the species with a lambda.x value over 10
dim(data_for_plot[data_for_plot$lambda.x > 10, ]) # 704 

data_for_plot[data_for_plot$lambda.x < 10, ]$lambda.x



# # Combine the plots
# combined_plot <- ggplot() +
#   # World map as background
#   annotation_custom(ggplotGrob(map_plot), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
#   # Overlay violin plot
#   annotation_custom(ggplotGrob(violin_plot), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

# # Plot
# print(combined_plot)


# # Combine the plots
# combined_plot_box <- ggplot() +
#   # World map as background
#   annotation_custom(ggplotGrob(map_plot), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
#   # Overlay violin plot
#   annotation_custom(ggplotGrob(boxplot), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

# # Plot
# print(combined_plot_box)




# # Combine the plots
# combined_plot_box <- ggdraw() +
#   draw_plot(map_plot, 0, 0, 1, 1) +
#   draw_plot(boxplot, 0, 0.19, 1, 0.51, hjust = 0, vjust = 0)

# # Plot
# print(combined_plot_box)

# # Combine the plots
# combined_plot_box <- ggdraw() +
#   draw_plot(map_plot, 0, 0, 1, 1) +
#   draw_plot(boxplot, 0, 0, 1, 0.9, hjust = 0, vjust = 0)

# # Plot
# print(combined_plot_box)


# Can I find all the species which are deemed as tropical but which are located north and south of the tropics which would be north of 23.5 and south of -23.5
tropical_species <- data_for_plot[data_for_plot$biome_affiliation == "TRF", ]
head(tropical_species)

# Just find the species which dont have a lat_band = -20_-10, -10_0, 0_10, 10_20
tropical_species_not_in_tropics <- tropical_species[!tropical_species$lat_band %in% c("-30_-20","-20_-10", "-10_0", "0_10", "10_20","20_30"), ]
head(tropical_species_not_in_tropics)
dim(tropical_species_not_in_tropics) # 1221 species. 
tropical_species_not_in_tropics$wcvp_taxon_name

# Can i now split these species into different latitude bands?
# Create a dataframe to store the data
data_for_plot_tropical_species_not_in_tropics <- data_for_plot[data_for_plot$wcvp_taxon_name %in% tropical_species_not_in_tropics$wcvp_taxon_name, ] # This holds all the rows where the species name is one of the problematic species names
head(data_for_plot_tropical_species_not_in_tropics)

# Now I want to find the rows in this dataset where the lat_band is not -30_-20, -20_-10, -10_0, 0_10, 10_20, 20_30
problems <- (data_for_plot_tropical_species_not_in_tropics[!data_for_plot_tropical_species_not_in_tropics$lat_band %in% c("-30_-20","-20_-10", "-10_0", "0_10", "10_20","20_30"),])
problems[,c(1,2)]
dim(problems)

data_for_plot_no_problem <- data_for_plot

# doing it with for loop
for (i in seq_along(problems$wcvp_taxon_name)) {
	problem_sp <- problems$wcvp_taxon_name[i]
	lat_band <- problems$lat_band[i]
	data_for_plot_no_problem <- subset(data_for_plot_no_problem, !(wcvp_taxon_name == problem_sp & lat_band == lat_band))
}

# check if the rows are removed
problems_still_in_data_for_plot <- problems[problems$wcvp_taxon_name %in% data_for_plot_no_problem$wcvp_taxon_name & problems$lat_band %in% data_for_plot_no_problem$lat_band, ]
problems_still_in_data_for_plot[,c(1,2)]



# update data_for_plot_counts aswell.
data_for_plot_counts_90 <- data_for_plot_no_problem %>%
	group_by(lat_band, in_tropical_rainforest90) %>%
	summarize(count = n())

data_for_plot_counts_biome_affiliation <- data_for_plot_no_problem %>%
	group_by(lat_band, biome_affiliation) %>%
	summarize(count = n())

head(data_for_plot_no_problem)


# Create the boxplot
boxplot_no_problems <- ggplot(data_for_plot_no_problem, aes(x = log(lambda.x), y = lat_band, fill = as.factor(biome_affiliation))) +
  geom_boxplot(alpha = 1) +
  scale_fill_manual(values = c("TRF" = cols[2], "Non TRF" = cols[1], "Widespread" = cols[3]),
                    labels = c("TRF" = "Tropical Rainforest", "Non TRF" = "Outside Tropical Rainforest", "Widespread" = "Widespread")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove panel background
    axis.text.x = element_blank() # Remove the x axis text as it is not needed
  ) +
  labs(x = expression("log((lineages/myr)"), y = "") +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key = element_rect(fill = NA, color = NA),
    legend.key.size = unit(1.5, "lines"),
    panel.background = element_rect(fill = "transparent", color = NA),
    axis.text.y = element_blank()
  ) +
  geom_text(data = data_for_plot_counts_biome_affiliation, aes(label = count, x = 9, y = lat_band, color = as.factor(biome_affiliation)),
            position = position_dodge(width = 0.9), size = 2.2) + # Map color to biome_affiliation
  scale_colour_manual(values = c("TRF" = cols[2], "Non TRF" = cols[1], "Widespread" = cols[3]))
boxplot_no_problems


# Combine the plots
combined_plot_box <- ggdraw() +
  draw_plot(map_plot, 0, 0, 1, 1) +
  draw_plot(boxplot_no_problems, 0, 0.05, 1, 0.85, hjust = 0, vjust = 0)


# Setting output folder
output_folder <- "/home/au543206/GenomeDK/Trf_models/workflow/05_figures" # Local
output_folder <- "/home/owrisberg/Trf_models/workflow/05_figures" # Srun



# Save the plot
pdf(file = file.path(output_folder,"Clads_world_map_test.pdf"),
    width = 10,
    height = 6)

print(combined_plot_box)

dev.off()


# Can I use the data to create a vertical density plot to show the latitudinal diversity gradient of the species.
# I can use the data_for_plot_no_problem dataframe to create a density plot of the species in each latitude band.

data_for_plot_counts_biome_affiliation <- data_for_plot_counts_biome_affiliation %>%
  mutate(lat_band = factor(lat_band, levels = sort(unique(lat_band))))

density_plot <- ggplot(data_for_plot_counts_biome_affiliation, aes(x = lat_band, y = count, fill = biome_affiliation)) +
  geom_bar(stat = "identity", color = "black", alpha = 1) +
  theme_void() +
  theme(
	legend.position = "none",
	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
  ) +
  labs(
    title = "",
    x = "",
    y = "",
    fill = ""
  ) +
  scale_fill_manual(values = c("TRF" = cols[2], "Non TRF" = cols[1], "Widespread" = cols[3])) +
  coord_flip() +
  scale_y_reverse()


# Can we add the density plot to the the combined_plot_box on the right side of the boxplot using cowplot

combined_plot_box_dens <- plot_grid(boxplot_no_problems, density_plot, ncol = 2, rel_widths = c(1, 0.2), align = "h", axis = "bt")
combined_plot_box_dens

# Now I need to draw the world map underneath the box plots on combined_plot_box_dens
combined_plot_box_dens <- ggdraw() +
  draw_plot(map_plot, 0, 0, 0.79, 1) +
  draw_plot(combined_plot_box_dens, 0, 0.35, 1, 0.31, hjust = 0, vjust = 0)

pdf(file = file.path(output_folder,"box_dens.pdf"),
	width = 10,
	height = 12)

print(combined_plot_box_dens)

dev.off()

# Combine the plots
combined_plot_box <- ggdraw() +
  draw_plot(map_plot, 0, 0, 1, 1) +
  draw_plot(boxplot, 0, 0.19, 1, 0.51, hjust = 0, vjust = 0)


# Can we do a simple test to see if there is a correlation between the lambda values and the latitude bands
# I would do this by creating a new dataframe which only includes the lambda values and the latitude bands
head(data_for_plot)
data_for_correlation <- data_for_plot[,c(1,2,3)]

# Can we now calculate the mean lambda value for each latitude band
mean_lambda <- data_for_correlation %>%
  group_by(lat_band) %>%
  summarize(mean_lambda = mean(lambda.x))

mean_lambda

# Can we convert the lat band to abselute latitude values
mean_lambda$lat_band_abs <- (c(40,30,20,10,0,0,10,20,30,40,50,60))

mean_lambda

# Is there a linear relationship between the latitude and the mean lambda values
cor_test <- cor(mean_lambda$mean_lambda, mean_lambda$lat_band_abs)
summary(lm(mean_lambda$mean_lambda ~ mean_lambda$lat_band_abs))

# Can we plot the mean lambda values against the latitude bands
ggplot(mean_lambda, aes(x = lat_band_abs, y = mean_lambda)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Absolute Latitude", y = "Mean Tip Speciation Rate") +
  theme_minimal() +
  theme(
	panel.background = element_blank()  # Remove panel background
  )

# Save as png
ggsave(file = file.path(output_folder, "mean_lambda_vs_latitude.png"), width = 10, height = 6)

output_folder
