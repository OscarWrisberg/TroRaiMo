#########################################################################################################################
################################################## Loading packages  ####################################################
#########################################################################################################################

# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("data.table", "ape", "phytools", "geiger", "dplyr", "ggplot2", "viridis","hrbrthemes", "cowplot", "MetBrewer","tidyverse", "posterior") #

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
folder_path <- "/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output/Esse_done_results/"

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

file_Nymphaeales <- "/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output/Esse_output_Nymphaeales_hidden_states_0.1.log"
file_Pandanales <- "/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output/Esse_output_Pandanales_hidden_states_0.1.log"
file_Piperales <- "/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output/Esse_output_Piperales_hidden_states_0.1.log"
file_Arecales <- "/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output/Esse_output_Arecales_hidden_states_0.1_10_flush.log"
file_Zingiberales <- "/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output/Esse_output_Zingiberales_hidden_states_0.1.log"


lambda_A_0 <- coda::as.mcmc(esse_raw$lambda_A_0)
lambda_B_0 <- coda::as.mcmc(esse_raw$lambda_B_0)
lambda_W_0 <- coda::as.mcmc(esse_raw$lambda_W_0)
lambda_A_1 <- coda::as.mcmc(esse_raw$lambda_A_1)
lambda_B_1 <- coda::as.mcmc(esse_raw$lambda_B_1)
lambda_W_1 <- coda::as.mcmc(esse_raw$lambda_W_1)

head(esse_raw_mcmc)

# Plot all the lambda values in a single plot
plot1 <- traceplot(lambda_A_0)
plot2 <- traceplot(lambda_B_0)
plot3 <- traceplot(lambda_W_0)
plot4 <- traceplot(lambda_A_1)
plot5 <- traceplot(lambda_B_1)
plot6 <- traceplot(lambda_W_1)

# Combine the plots
plot_grid(plot1, plot2, plot3, plot4, plot5, plot6)

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


#################################################################
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Sample data based on your description
data <- data.frame(
  order = c(rep("Alismatales", 6), rep("Poales", 6)),
  Var = rep("lambda", 12),
  State = rep(c("A", "B", "W"), 4),
  Hiddenstate = rep(c(0, 0, 0, 1, 1, 1), 2),
  value = c(0.364247, 0.246137, 0.009381, 0.389014, 0.533666, 0.000066,
            0.274837, 0.346137, 0.013481, 0.289014, 0.433666, 0.000086)
)

# Plot using ggplot2
ggplot(Esse_dataframe, aes(x = value, y = State, fill = factor(Hiddenstate))) +
  geom_violin(trim = FALSE) +
#  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  labs(title = "Violin plot of Lambda values",
       x = "Lambda Variables",
       y = "Estimated Value",
       fill = "Hiddenstate") +
	   xlim(0, 20) +
  theme_minimal()

Esse_nymp_low_test <- Esse_dataframe[which(Esse_dataframe$order == "Nymphaeales"),]
# removing the values that are over 50
Esse_nymp_low_test <- Esse_nymp_low_test[which(Esse_nymp_low_test$value < 5),]


sorted_vector <- sort(Esse_dataframe[which(Esse_dataframe$order == "Nymphaeales"),"value"])
sorted_vector


# Plot using ggplot2
ggplot(Esse_dataframe, aes(x = value, y = State, fill = factor(Hiddenstate))) +
  geom_violin(trim = TRUE) +
  #geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  labs(title = "Violin plot of Lambda values",
       x = "Lambda Variables",
       y = "Estimated Value",
       fill = "Hiddenstate") +
	xlim(0, 5) +
  theme_minimal()

install.packages("ggridges")
library(ggridges)

ggplot(Esse_dataframe, aes(x = value, y = State, fill = State)) +
  geom_density_ridges(alpha = 0.8, scale = 1) +
  labs(title = "Density Ridge plot of Lambda values",
       x = "Lambda Variables",
       y = "Proportion",
       fill = "State") +
  xlim(0, 20) +
  theme_minimal()



ggplot(Esse_dataframe, aes(x = value, fill = State)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density plot of Lambda values",
       x = "Lambda Variables",
       y = "Density",
       fill = "State") +
  xlim(0, 20) +
  theme_minimal()

Esse_dataframe_0 <- Esse_dataframe[which(Esse_dataframe$Hiddenstate == 0),]
Esse_dataframe_1 <- Esse_dataframe[which(Esse_dataframe$Hiddenstate == 1),]


################################

# Plot using ggplot2
ggplot() +
	# Top density plot
	geom_density(data = Esse_dataframe_0, aes(x = value, fill = State), alpha = 0.5) +
	# Bottom density plot (mirrored)
	geom_density(data = Esse_dataframe_1, aes(x = value, y = ..density.. * -1, fill = State), alpha = 0.5) +
	labs(title = "Mirrored Density plot of Lambda values",
			 x = "Lambda Variables",
			 y = "") +
	xlim(0, 20) +
	theme_minimal()


ggplot() +
  # Top density plot
  geom_density(data = Esse_dataframe_0, aes(x = value, fill = State), alpha = 0.4) +
  # Bottom density plot (mirrored)
  geom_density(data = Esse_dataframe_1, aes(x = value, y = -..density.., fill = State), alpha = 0.4) +
  scale_y_continuous(
    name = "Density",
    sec.axis = sec_axis(~., name = "Density")
  ) +
  labs(title = "Mirrored Density plot of Lambda values",
       x = "Lambda Variables",
       y = "Density",
       fill = "State") +
  xlim(0,100) +
  theme_minimal() +
  theme(
    axis.title.y.left = element_blank(),
    axis.text.y.left = element_blank(),) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid")

# I want to so