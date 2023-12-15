#########################################################################################################################
################################################## Loading packages  ####################################################
#########################################################################################################################

# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("raster","ncdf4","ggplot2","dplyr","data.table")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


#########################################################################################################################
################################################-- Loading Data --#######################################################
#########################################################################################################################

# # Command line arguments
data_dir <- commandArgs(trailingOnly = TRUE)[1]
output <- commandArgs(trailingOnly = TRUE)[2]


#paleo_clim_folder <- "/home/au543206/GenomeDK/Trf_models/data/paleo_clim"
paleo_clim_folder <- data_dir

# Here I list the order of the paleo_clim files
paleo_clim_order <- c(
"540Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"520Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"500Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"480Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"460Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"440Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"420Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"400Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"380Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"360Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"340Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"320Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"300Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"280Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"260Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"240Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"220Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"200Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"180Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"160Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"140Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"120Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"100Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"80Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"60Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"40Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"20Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
"0Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv"
)

#Constructing a data frame to store the data
time_series_data <- data.frame(time_ago = numeric(),
	total_area_rainforest = numeric(),
	total_area_other = numeric())

# I am using this for loop to go through each of the raster files.
# For each of the raster files I want to extract the area for the cells thar are defined as Tropical Rainforest (12 and 13) and the area for the other cells.
# All these variables should be saved to the time_series_data data frame with one column being the time and the other columns being the variables
# The time is gathered from the file name and is the first element of the file name when it is split by "_"

for (i in seq_along(paleo_clim_order)){
	# Read the raster file
	paleo_file <- fread(file.path(paleo_clim_folder, paleo_clim_order[i]), header=TRUE, sep=",",verbose=FALSE, quote="")

	# Extract the time from the file name
	time <- strsplit(paleo_clim_order[i], "_")[[1]][1]
	print(time)
	# Remove the "Ma" from the time
	time <- gsub("Ma", "", time)

	# Calculate the total area of the cells that are defined as tropical and others
	total_area_rainforest <- sum(paleo_file[paleo_file$koppen == 12 | paleo_file$koppen == 13]$area)
	total_area_other <- sum(paleo_file[paleo_file$koppen != 12 | paleo_file$koppen != 13]$area)
	
	# I think this is fine for now
	# Add the vars to the time_series_data data frame
	time_series_data <- rbind(time_series_data, data.frame(time_ago = time, total_area_rainforest = total_area_rainforest, total_area_other = total_area_other))

}


#time_series_data

# Convert time_ago to numeric
time_series_data$time_ago <- as.numeric(time_series_data$time_ago)

# Plot the data with ggplot Tropical area should be green dots and the other area should be ivory dots
# ggplot(data = time_series_data) +
#   geom_point(aes(x = time_ago, y = total_area_rainforest), color = "green") +
#   geom_point(aes(x = time_ago, y = total_area_other), color = "black") +
#   geom_smooth(aes(x = time_ago, y = total_area_rainforest), color = "green") +
#   geom_smooth(aes(x = time_ago, y = total_area_other), color = "black") +
#   theme_bw()

write.table(time_series_data, file = output, sep = ",", row.names = FALSE)

