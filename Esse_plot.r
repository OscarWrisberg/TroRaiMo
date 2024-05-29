#########################################################################################################################
################################################## Loading packages  ####################################################
#########################################################################################################################

# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("data.table", "ape", "phytools", "geiger", "dplyr", "ggplot2", "viridis","hrbrthemes", "cowplot", "MetBrewer","tidyverse", "coda") #

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

#########################################################################################################################
################################################## Loading Esse data ####################################################
#########################################################################################################################

# Loading data path
folder_path <- "/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output"

# Get a list of all Rdata files in the folder
file_list <- list.files(folder_path, pattern = "Esse_output_.*\\_hidden_states_0.1.log", full.names = TRUE)


# Create dataframe for tips and tip rate speciation
Esse_dataframe <- data.frame()

  
# Loop through each file in the folder
for (i in 1:length(file_list)) { 

	# Keeping track of the progress
	print(paste("Processing file", i, "of", length(file_list)))

	# Extractinc the file name
	file = file_list[i]

	# Extract the order name from the file name
	order_name <- gsub("Esse_output_(.*)\\_hidden_states_0.1.log","\\1", basename(file))
	
	print(order_name)
	# Load the Rdata file
	esse_raw <- fread(file)

	# Convert the wide data to long format
	long_format_ESSE <- esse_raw %>%
	  pivot_longer(cols = everything(), names_to = c("Var", "State", "Hiddenstate"), names_sep = "_")

	# Create a subset which contains only the lambda values
	lambda_mu_values <- long_format_ESSE %>%
	  filter(Var == "lambda" | Var == "mu")

	print(lambda_mu_values)

	# Creating a dataframe with the order name
	order_data <- data.frame(order = order_name,Var = lambda_mu_values$Var,State = lambda_mu_values$State,Hiddenstate = lambda_mu_values$Hiddenstate, value = lambda_mu_values$value)

	# Append the tip names and tip rate speciation to the dataframe
	Esse_dataframe <- rbind(Esse_dataframe, order_data)
	
	# print dim of the dataframe to keep track of the progress
	#cat("The dataset contains ",dim(lambda_mu_values)[1], " rows and ", dim(lambda_mu_values)[2], " columns\n")

	# Remove the esse_raw object from the environment
	rm(esse_raw)
}

head(Esse_dataframe)

esse_raw <- fread(file)
esse_raw
esse_raw_mcmc <- coda::as.mcmc(esse_raw)
ess <- effectiveSize(esse_raw_mcmc)
ess
#########################################################################################################################
################################################## Plotting Esse data ####################################################
#########################################################################################################################



# Here I want to plot a boxplot of the lambda_B_1 values

trfcol <- "chartreuse"
non_trfcol <- "deeppink"
widespread <- "aquamarine2"

ggplot(lambda_values, aes(x = State, y = value)) +
	geom_boxplot(fill = c(non_trfcol, trfcol,widespread), color = "black") +
	scale_x_discrete(labels = c("Not in TRF", "In TRF")) +
	ylab("Lambda") +
	xlab("") +
	theme_ipsum() + 
	theme(legend.position = "none")

ggplot(lambda_values, aes(x = value, fill = State)) +
	geom_density(alpha = 0.5) +
	scale_fill_manual(values = c(non_trfcol, trfcol, widespread), labels = c("Not in TRF", "In TRF", "Widespread")) +
	ylab("Density") +
	xlab("Lambda") +
	theme_ipsum() +
	labs(title = "Lambda") +
	scale_x_log10()

ggplot(lambda_values, aes(x = value, fill = State)) +
	geom_density(alpha = 0.5) +
	scale_fill_manual(values = c(non_trfcol, trfcol, widespread), labels = c("Not in TRF", "In TRF", "Widespread")) +
	ylab("Density") +
	xlab("Lambda") +
	theme_ipsum() +
	labs(title = "Mu") +
	scale_x_log10()
