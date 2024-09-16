#########################################################################################################################
################################################## Loading packages  ####################################################
#########################################################################################################################

# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("data.table", "ape", "phytools", "geiger", "dplyr", "ggplot2", "viridis", "cowplot", "MetBrewer","caper","broom", "phylolm") #

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Setting output folder
output_folder <- "/home/au543206/GenomeDK/Trf_models/workflow/05_figures" # local
output_folder <- "/home/owrisberg/Trf_models/workflow/05_figures" # remote

########################################################################################################################
############################################## Checking for duplicates  ################################################
########################################################################################################################

# Clads on orders
folderpath_1 <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/orders" # local
folderpath_1 <- "/home/owrisberg/Trf_models/workflow/02_adding_orders/pruning/orders" # remote

file_list_1 <- list.files(folderpath_1, pattern = "Clads_output_.*\\.Rdata", full.names = TRUE)

#Clads on suborders
folderpath_2 <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/subset_of_orders" # local
folderpath_2 <- "/home/owrisberg/Trf_models/workflow/02_adding_orders/pruning/subset_of_orders" # remote
file_list_2 <- list.files(folderpath_2, pattern = "Clads_output_.*\\.Rdata", full.names = TRUE)

# join the file lists
file_list <- c(file_list_1, file_list_2)

# Length of file list
length(file_list) # 291

#########################################################################################################################
############################################ Loading ClaDs runs on orders  ##############################################
#########################################################################################################################

# Folder path for distribution data
folder_path <- "/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data" # local
folder_path <- "/home/owrisberg/Trf_models/workflow/03_distribution_data" # remote

# Get a list of all distribution files in the folder
distribution_files <- list.files(folder_path, pattern = ".*_distribution_data_ClaDs.txt", full.names = TRUE)
distribution_files

# Correct column names
correct_names <- c("wcvp_taxon_name", "occurrences_trf","occurrences_non_trf","proportion_in_tropical_rainforest","proportion_outside_tropical_rainforest")

# Set the theshold values
thresholds <- c(0.33)

# Create a list of comparative.dataframes
cdat_list <- list()

# Files where not all species are in ditribution data
broken_runs_dist <- c()

# list of broken ClaDs runs
broken_runs <- c()

# I need to convert this part to a for loop which loops through each Clads file in the folder and create a comparative.dataframe for each ClaDs run.
for (i in seq_along(file_list)){

	# Keeping track of the progress
	print(paste("Processing file", i, "of", length(file_list)))

	file <- file_list[i]

	print(paste0("Order: ",gsub("Clads_output_(.*)\\.Rdata", "\\1", basename(file))))

	load(file)

	# Selecting the tree and renaming the tips
	tree <- CladsOutput$tree
	tree$tip.label <- gsub("_", " ", tree$tip.label)

	# Selecting the tip rate speciation
	CladsOutput$lambdatip_map

	# Selecting the tips and renaming them
	taxa <- tree$tip.label
	taxa <- gsub("_", " ", taxa)

	# I also need to find the biome states of the tips in the tree..
	order <- gsub("Clads_output_(.*)\\.Rdata", "\\1", basename(file))
	order

	# Select the distribution file for the order
	distribution_file <- distribution_files[grep(paste0(order,"_"), distribution_files)]
	distribution_file

	if(length(distribution_file) > 1) {
		print("Error: More than one distribution file found")
		print(distribution_file)
		break
	}

	if(length(distribution_file) == 0) {
		broken_runs <- c(broken_runs, order)
		print("Error: No distribution file found")
		next
	}

	# Load the distribution data
	distribution_data <- fread(distribution_file)


	if(all(names(distribution_data) %in% correct_names)) {
		# Do nothing	
	} else {
		# Rename the columns
		setnames(distribution_data, c("wcvp_taxon_name", "occurrences_trf","occurrences_non_trf","proportion_in_tropical_rainforest","proportion_outside_tropical_rainforest"))
	}

	# Iterate over each threshold value to create the corresponding columns
	for (threshold in thresholds) {
	  # Create the column names based on the threshold
	  in_trf_col <- paste0("in_tropical_rainforest", threshold * 100)
	  outside_trf_col <- paste0("outside_tropical_rainforest", threshold * 100)
	  biome_affiliation_col <- paste0("biome_affiliation_", threshold * 100)
	
	  # Add columns to the distribution data
	  distribution_data[[in_trf_col]] <- distribution_data$proportion_in_tropical_rainforest > threshold
	  distribution_data[[outside_trf_col]] <- distribution_data$proportion_outside_tropical_rainforest > threshold
	
	  # Add the biome affiliation column based on the conditions
	  distribution_data[[biome_affiliation_col]] <- ifelse(
	    distribution_data[[in_trf_col]] == TRUE & distribution_data[[outside_trf_col]] == FALSE,
	    "TRF",
	    ifelse(
	      distribution_data[[outside_trf_col]] == TRUE & distribution_data[[in_trf_col]] == FALSE,
	      "Non TRF",
	      ifelse(
	        distribution_data[[in_trf_col]] == TRUE & distribution_data[[outside_trf_col]] == TRUE,
	        "Widespread",
	        "Error"
	      )
	    )
	  )
	}

	# Create a dataframe with the tip labels and file names
	dat <- data.frame(
	taxa = taxa,
	lambda = CladsOutput$lambdatip_map,
	extinction = CladsOutput$eps_map,
	order = order
	)

	# Check if all the tips in dat are in the distribution data
	if(!all(dat$taxa %in% distribution_data$wcvp_taxon_name)) {
		print("Error: Not all tips in the tree are in the distribution data")
		# Dropping the taxa which are not in the distribution data
		dat <- dat[dat$taxa %in% distribution_data$wcvp_taxon_name,]
		# Dropping distribution data taxa which are not in dat
		distribution_data <- distribution_data[distribution_data$wcvp_taxon_name %in% dat$taxa,]
	}

	if(!all(distribution_data$wcvp_taxon_name %in% dat$taxa)) {
		print("Error: Not all tips in the distribution data are in the tree")
		# Dropping the taxa which are not in the tree
		distribution_data <- distribution_data[distribution_data$wcvp_taxon_name %in% dat$taxa,]
		# Dropping dat taxa which are not in distribution data
		dat <- dat[dat$taxa %in% distribution_data$wcvp_taxon_name,]
	}

	if(nrow(dat) == 0) {
		print("Error: Empty dataframe")
		next
	}

	# add the biome affiliation to the dataframe
	dat <- cbind(dat, distribution_data[match(dat$taxa, distribution_data$wcvp_taxon_name), c("biome_affiliation_33")])

	# Create a comparative.data object
	cdat <- comparative.data(data = dat, phy = tree, names.col = "taxa")

	print(dim(cdat$data))

	# Append the cdat object to the list
	cdat_list[[i]] <- cdat
}

# How many of the elements in the cdat_list are not null
sum(!sapply(cdat_list, is.null)) # 263

# save the cdat_list
save(cdat_list, file = file.path(output_folder, "cdat_list.rda"))
load(file.path(output_folder, "cdat_list.rda")) 

##########################################################################################################################
###################################################---- Fitting the PGLS ----#############################################
##########################################################################################################################

# Testing the pgls from caper
significant_orders <- list()
non_significant_orders <- list()
significant_models <- list()
non_significant_models <- list()
only_one_biome_affiliation <- list()
broken_runs <- list()

length(cdat_list)

for (i in seq_along(cdat_list)){
	cdat <- cdat_list[[i]]

	# keep track of progress
	print(paste("Processing file", i, "of", length(cdat_list)))
	order <- unique(cdat$data$order)

	if(is.null(cdat)) {
		print("Null")
		broken_runs <- append(broken_runs, order)
		next
	}

	if(length(cdat$phy$tip.label) == 0) {
		print("Empty dataframe")
		broken_runs <- append(broken_runs, order)
		next
	}

	if(length(unique(cdat$data$biome_affiliation_33)) == 1) {
		print("Only one biome affiliation")
		only_one_biome_affiliation <- c(only_one_biome_affiliation, order)
		next
	}

	output_test <- pgls(lambda ~ biome_affiliation_33, data = cdat)

	output_summary <- list(summary(output_test))


	 if( any(output_summary[[1]]$coefficients[,4] < 0.05) ) {
	 	print("Significant")
		significant_orders <- append(significant_orders, order)
		significant_models <- append(significant_models, output_summary)
	 } else {
	 	print("Not significant")
		non_significant_orders <- append(non_significant_orders, order)
		non_significant_models <- append(non_significant_models, output_summary)
	 }
}

# Now we need to save the different outputs so we can work on them locally.

# Save the significant orders
save(significant_orders, file = file.path(output_folder, "significant_orders_pgls.rda"))

# Save the non-significant orders
save(non_significant_orders, file = file.path(output_folder, "non_significant_orders_pgls.rda"))

# Save the significant models
save(significant_models, file = file.path(output_folder, "significant_models_pgls.rda"))

# Save the non-significant models
save(non_significant_models, file = file.path(output_folder, "non_significant_models_pgls.rda"))

# Save the broken runs
save(broken_runs, file = file.path(output_folder, "broken_runs_pgls.rda"))

# Save the only one biome affiliation
save(only_one_biome_affiliation, file = file.path(output_folder, "only_one_biome_affiliation_pgls.rda"))

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

# Load the different outputs
load(file.path(output_folder, "significant_orders_pgls.rda"))
load(file.path(output_folder, "non_significant_orders_pgls.rda"))
load(file.path(output_folder, "significant_models_pgls.rda"))
load(file.path(output_folder, "non_significant_models_pgls.rda"))
load(file.path(output_folder, "broken_runs_pgls.rda"))
load(file.path(output_folder, "only_one_biome_affiliation_pgls.rda"))

# Print the results
print(significant_orders)

print(non_significant_orders)

print(significant_models)

print(non_significant_models)

print(only_one_biome_affiliation)

print(broken_runs)



# Can i loop through all the significant models and  just get the coefficient table?
coefficient_table_significant <- lapply(signifcant_models, function(x) x$coefficients)

significant_models[[1]]$coefficients

# Can we loop throigh all the significant models and find out which of the variables are significant?
significant_variables <- lapply(significant_models, function(x) x$coefficients[,4] < 0.05)
significant_variables

# Initialize a named vector to store the counts
significance_counts <- c(
  "(Intercept)" = 0,
  "biome_affiliation_33TRF" = 0,
  "biome_affiliation_33Widespread" = 0
)

# Iterate through each list of significant variables and count TRUE occurrences
for (model in significant_variables) {
  for (variable in names(model)) {
    if (model[[variable]] == TRUE) {
		cat("Variable", variable, "is significant\n")
    	significance_counts[variable] <- significance_counts[variable] + 1
    }
  }
}

significance_counts # is it really 14/14/14


# I need to find out how many of these have a positive and negative coefficient
significance_counts_positive <- c(
  "(Intercept)" = 0,
  "biome_affiliation_33TRF" = 0,
  "biome_affiliation_33Widespread" = 0
)

significance_counts_negative <- c(
  "(Intercept)" = 0,
  "biome_affiliation_33TRF" = 0,
  "biome_affiliation_33Widespread" = 0
)

# Iterate through each list of significant variables and count TRUE occurrences
for (model in significant_models) {
  for (variable in names(model$coefficients[,4])) {
	if (model$coefficients[variable,4] < 0.05) {
	  if (model$coefficients[variable,1] > 0) {
		significance_counts_positive[variable] <- significance_counts_positive[variable] + 1
	  } else {
		significance_counts_negative[variable] <- significance_counts_negative[variable] + 1
	  }
	}
  }
}

significance_counts_positive

significance_counts_negative

#########################################################################################################################
###################################################---- Pgls on entire tree ----#############################################
#########################################################################################################################

# Can I instead combine all the CladS Outputs and then prune the aligned Smith and Brown tree to see if I get a different result?
distribution_data_merged <- readRDS("/home/au543206/GenomeDK/Trf_models/workflow/04_results/distribution_data_merged.rds") # This dataset is created by the Clads_plot.r script
distribution_data_merged <- readRDS("/home/owrisberg/Trf_models/workflow/04_results/distribution_data_merged.rds") # This dataset is created by the Clads_plot.r script

# Now I need to load the smith and brown tree
smb_tree <- read.tree("/home/au543206/GenomeDK/Trf_models/data/GBMB_pruned.tre")
smb_tree <- read.tree("/home/owrisberg/Trf_models/data/GBMB_pruned.tre")

# change _ to " " in the tip labels
smb_tree$tip.label <- gsub("_", " ", smb_tree$tip.label)

# Are there any duplicates in the smb tree?

# How many of the tip in the distribution data are in the tree?
sum(distribution_data_merged$tip_label %in% smb_tree$tip.label) # 0

# Remove the species which are not in the distribution data
smb_tree_pruned <- drop.tip(smb_tree, setdiff(smb_tree$tip.label, distribution_data_merged$tip_label))

# are there any duplicates in the distribution data?
sum(duplicated(distribution_data_merged$tip_label)) # 0
length(distribution_data_merged[which(duplicated(distribution_data_merged$tip_label)),]$tip_label)


any(duplicated(distribution_data_merged$tip_label)) # FALSE
any(duplicated(smb_tree_pruned$tip.label)) # FALSE

smb_tree_pruned$node.label <- NULL

# Create a comparative.data object
cdat_smb <- comparative.data(data = distribution_data_merged, phy = smb_tree_pruned, names.col = "tip_label")

dim(cdat_smb$data)
head(cdat_smb$data)

########################################################################################################
########---- OBS these models takes a very long time to run and more than 100G of RAM ----##############
########################################################################################################

# Fit the PGLS model
# output_test_smb_10 <- pgls(lambda ~ biome_affiliation_10, data = cdat_smb)
# output_summary_smb_10 <- summary(output_test_smb_10)

# output_test_smb_25 <- pgls(lambda ~ biome_affiliation_25, data = cdat_smb)
# output_summary_smb_25 <- summary(output_test_smb_25)

# output_test_smb_33 <- pgls(lambda ~ biome_affiliation_33, data = cdat_smb)
# output_summary_smb_33 <- summary(output_test_smb_33)

# output_test_smb_50 <- pgls(lambda ~ biome_affiliation_50, data = cdat_smb)
# output_summary_smb_50 <- summary(output_test_smb_50)

########################################################################################################
#############################################---- PhyloLM ----##########################################
########################################################################################################


head(cdat_smb)

smb_tree_pruned$node.label <- NULL
head(distribution_data_merged)

options(future.globals.maxSize= 1048576000)

# Changing the row names to be the tip_labels
rownames(distribution_data_merged) <- distribution_data_merged$tip_label

# Fit the phylolm model
output_test_smb_phylolm_10 <- phylolm(lambda ~ biome_affiliation_10, data = distribution_data_merged, phy = smb_tree_pruned, boot = 1000)
summary(output_test_smb_phylolm_10)

output_test_smb_phylolm_25 <- phylolm(lambda ~ biome_affiliation_25, data = distribution_data_merged, phy = smb_tree_pruned, boot = 1000)
summary(output_test_smb_phylolm_25)

output_test_smb_phylolm_33 <- phylolm(lambda ~ biome_affiliation_33, data = distribution_data_merged, phy = smb_tree_pruned, boot = 1000)
summary(output_test_smb_phylolm_33)

output_test_smb_phylolm_40 <- phylolm(lambda ~ biome_affiliation_40, data = distribution_data_merged, phy = smb_tree_pruned, boot = 1000)
summary(output_test_smb_phylolm_40)

class(summary(output_test_smb_phylolm_10)[2])

# Function to convert scientific notation to decimal format
convert_to_decimal <- function(x) {
  if (is.numeric(x)) {
    format(x, scientific = FALSE)
  } else {
    x
  }
}

# Function to recursively apply conversion to lists and data frames
convert_summary <- function(x) {
  if (is.list(x)) {
    lapply(x, convert_summary)
  } else if (is.data.frame(x)) {
    as.data.frame(lapply(x, function(col) {
      if (is.numeric(col)) {
        format(col, scientific = FALSE)
      } else {
        col
      }
    }))
  } else {
    convert_to_decimal(x)
  }
}

# Apply conversion to the summary object
summary_converted_10 <- convert_summary(summary(output_test_smb_phylolm_10)[2])
summary_converted_25 <- convert_summary(summary(output_test_smb_phylolm_25)[2])
summary_converted_33 <- convert_summary(summary(output_test_smb_phylolm_33)[2])
summary_converted_40 <- convert_summary(summary(output_test_smb_phylolm_40)[2])


# Print or inspect the converted summary
summary_converted_10
summary_converted_25
summary_converted_33
summary_converted_40

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

# Here I also want to test if there is a difference in sampling rates for the phylogenies which are significant and non-significant
# The first thing I will need to do is to find the sampling rate of all the trees which have a significant result.
significant_orders

# Find the sampling fractions for all these orders
list_of_sampling_fractions <- list.files("/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data", pattern = ".*_sampling_fraction\\.txt", full.names = TRUE)
list_of_sampling_fractions <- list.files("/home/owrisberg/Trf_models/workflow/03_distribution_data", pattern = ".*_sampling_fraction\\.txt", full.names = TRUE)

# Remove the sampling fractions which have biome in them
list_of_sampling_fractions <- list_of_sampling_fractions[!grepl("biome", list_of_sampling_fractions)]

# Remove all the sampling fractions which have sub_phylo in them
list_of_sampling_fractions <- list_of_sampling_fractions[!grepl("sub_phylo", list_of_sampling_fractions)]

# sampling fractions of the significant orders
sampling_fractions_significant <- list_of_sampling_fractions[grep(paste(paste0(significant_orders,"_"), collapse = "|"), list_of_sampling_fractions)]

sampling_fractions_significant

# Non significant orders
non_significant_orders

# Find the sampling fractions for all these orders
sampling_fractions_not_significant <- list_of_sampling_fractions[grepl(paste(paste0(non_significant_orders,"_"), collapse = "|"), list_of_sampling_fractions)]

# So now we have a list of all the sampling fractions for the significant and non-significant orders. 
#Now we need to load these files

# Load the sampling fractions
sampling_fractions_significant_data <- lapply(sampling_fractions_significant, fread)
sampling_fractions_significant_data


# Load the sampling fractions
sampling_fractions_not_significant_data <- lapply(sampling_fractions_not_significant, fread)
sampling_fractions_not_significant_data

# So what I need to do is calculate the sampling fraction for each of those lists.
# I guess what I have to do is select each of the genera in each list then select their sampling fraction and the number of species in that genus and that would give me the
# the total number of species in that genus

# This is a list of all the species in the significant orders
unlist(lapply(sampling_fractions_significant_data, function(x) x$species))

# Now I split the species names and take the first element
calculate_phylogeny_sampling_fraction <- function(df) {
  # Group by genus and summarize to get the effective number of sampled species
  df <- df %>%
    mutate(Genus = sapply(strsplit(as.character(species), " "), `[`, 1)) %>%
    group_by(Genus) %>%
    summarize(Species_Count = n(),
              Sampling_Fraction = first(sampling_freq)) %>%
    mutate(Effective_Sampled = Species_Count * Sampling_Fraction)
  
  # Calculate the overall sampling fraction
  total_effective_sampled <- sum(df$Effective_Sampled)
  total_species <- sum(df$Species_Count)
  
  overall_sampling_fraction <- total_effective_sampled / total_species
  
  return(overall_sampling_fraction)
}

# Assuming `phylogeny_list` is a list of dataframes
overall_sampling_fractions_not_significant <- lapply(sampling_fractions_not_significant_data, calculate_phylogeny_sampling_fraction)
overall_sampling_fractions_not_significant <- unlist(overall_sampling_fractions_not_significant)
overall_sampling_fractions_not_significant

# Assuming `phylogeny_list` is a list of dataframes
overall_sampling_fractions_significant <- lapply(sampling_fractions_significant_data, calculate_phylogeny_sampling_fraction)
overall_sampling_fractions_significant <- unlist(overall_sampling_fractions_significant)
overall_sampling_fractions_significant

shapiro.test(overall_sampling_fractions_significant)
shapiro.test(overall_sampling_fractions_not_significant)

t.test(overall_sampling_fractions_significant, overall_sampling_fractions_not_significant, alternative = "two.sided")

boxplot(overall_sampling_fractions_significant, overall_sampling_fractions_not_significant, names = c("Significant", "Non-significant"), ylab = "Sampling Fraction", xlab = "Significance")

# Can i make a ggplot violin plot instead?
sampling_fractions <- data.frame(
  Sampling_Fraction = c(overall_sampling_fractions_significant, overall_sampling_fractions_not_significant),
  Significance = rep(c("Significant", "Non-significant"), c(length(overall_sampling_fractions_significant), length(overall_sampling_fractions_not_significant)))
)

ggplot(sampling_fractions, aes(x = Significance, y = Sampling_Fraction, fill = Significance)) +
  geom_violin(trim = FALSE) +
  labs(title = "Sampling Fraction of Phylogenies",
	   x = "Significance",
	   y = "Sampling Fraction") +
  theme_minimal()


#################################

# Selecting colours for the plot
cols <- met.brewer("Morgenstern", n=7, type = "discrete")
cols
cols <- cols[c(3, 6)]

# Create the data frame
sampling_fractions <- data.frame(
  Sampling_Fraction = c(overall_sampling_fractions_significant, overall_sampling_fractions_not_significant),
  Significance = rep(c("Significant", "Non-significant"), c(length(overall_sampling_fractions_significant), length(overall_sampling_fractions_not_significant)))
)

# Plot with boxplot and jitter
ggplot(sampling_fractions, aes(x = Significance, y = Sampling_Fraction, fill = Significance)) +
  geom_boxplot(alpha = 0.9, fill = cols) +
  geom_jitter(color = "black", size = 1, alpha = 0.6, width = 0.2) +
  labs(
       x = "Significance",
       y = "Sampling Fraction") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11)
  )
# save the plot
plot_folder <- "/home/au543206/GenomeDK/Trf_models/workflow/05_figures"
plot_folder <- "/home/owrisberg/Trf_models/workflow/05_figures"
ggsave(file.path(plot_folder, "sampling_fraction_plot.png"), width = 6, height = 4, dpi = 300)

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
# Can we now see if the lower sampling fraction is simply due to the the size of the phylogenies?
# We need to count the number of tips per genera in the phylogenies so we know how many species are in each tree
# We can then combine the sampling fractions per genera and the number of species per genera to see if the true number of species in each tree is also different between the significant and non-significant trees


# Find all unique genera in the significant orders
# Extract species from each data frame
all_species_significant <- unlist(lapply(sampling_fractions_significant_data, function(x) x$species))

# Split species names into parts and extract the genus (first part)
# Assuming genus and species are separated by a space
genus_significant <- sapply(strsplit(all_species_significant, " "), `[`, 1)

# Get unique genera
unique_genera_significant <- unique(genus_significant)
unique_genera_significant

# Find all unique genera in the non-significant orders
# Extract species from each data frame
all_species_not_significant <- unlist(lapply(sampling_fractions_not_significant_data, function(x) x$species))

# Split species names into parts and extract the genus (first part)
# Assuming genus and species are separated by a space
genus_not_significant <- sapply(strsplit(all_species_not_significant, " "), `[`, 1)

# Get unique genera
unique_genera_not_significant <- unique(genus_not_significant)
unique_genera_not_significant

# Now we need to count the number of species in each genus for the significant and non-significant orders
# Function to count the number of species in each genus
count_species_per_genus <- function(df) {
  # Split species names into parts and extract the genus (first part)
  # Assuming genus and species are separated by a space
  df$Genus <- sapply(strsplit(as.character(df$species), " "), `[`, 1)
  
  # Count the number of species in each genus
  species_per_genus <- df %>%
	group_by(Genus) %>%
	summarize(Species_Count = n())
  
  return(species_per_genus)
}

# Count the number of species in each genus for the significant orders	
species_per_genus_significant <- lapply(sampling_fractions_significant_data, count_species_per_genus)
species_per_genus_significant

# Count the number of species in each genus for the non-significant orders
species_per_genus_not_significant <- lapply(sampling_fractions_not_significant_data, count_species_per_genus)
species_per_genus_not_significant

# Now we need to combine the sampling fractions and the number of species per genus
# Function to combine sampling fractions and number of species per genus
combine_sampling_fractions_species_per_genus <- function(sampling_df, species_df) {
  # Merge the two data frames by genus
  combined_df <- merge(sampling_df, species_df, by = "Genus", all = TRUE)
  
  # Replace NA values with 0
  combined_df[is.na(combined_df)] <- 0
  
  return(combined_df)
}

# Convert the sampling_fractions_data 
# find the sampling fractions of all the unique genera in the significant orders
sampling_fractions_significant_data <- lapply(sampling_fractions_significant_data, function(x) {
  x$Genus <- sapply(strsplit(as.character(x$species), " "), `[`, 1)
  x
})

# Convert the sampling_fractions_data
sampling_fractions_not_significant_data <- lapply(sampling_fractions_not_significant_data, function(x) {
  x$Genus <- sapply(strsplit(as.character(x$species), " "), `[`, 1)
  x
})

# Create a dataframe which contains all the sampling fractions and the number of species per genus
combined_sampling_fractions_species_per_genus_significant <- Map(combine_sampling_fractions_species_per_genus, sampling_fractions_significant_data, species_per_genus_significant)
combined_sampling_fractions_species_per_genus_not_significant <- Map(combine_sampling_fractions_species_per_genus, sampling_fractions_not_significant_data, species_per_genus_not_significant)

# Drop the species column
combined_sampling_fractions_species_per_genus_significant <- lapply(combined_sampling_fractions_species_per_genus_significant, function(x) x[,-2])
combined_sampling_fractions_species_per_genus_not_significant <- lapply(combined_sampling_fractions_species_per_genus_not_significant, function(x) x[, -2])

# Remove duplicate rows
combined_sampling_fractions_species_per_genus_significant <- lapply(combined_sampling_fractions_species_per_genus_significant, function(x) unique(x))
combined_sampling_fractions_species_per_genus_not_significant <- lapply(combined_sampling_fractions_species_per_genus_not_significant, function(x) unique(x))

# Now we add a third column which is the sampling_freq * species_count = total_species
combined_sampling_fractions_species_per_genus_significant <- lapply(combined_sampling_fractions_species_per_genus_significant, function(x) {
  x$Total_Species <- x$Species_Count / x$sampling_freq
  x
})

combined_sampling_fractions_species_per_genus_not_significant <- lapply(combined_sampling_fractions_species_per_genus_not_significant, function(x) {
  x$Total_Species <-  x$Species_Count / x$sampling_freq
  x
})

combined_sampling_fractions_species_per_genus_significant

# So now I calculate the total species count per dataframe in the list
total_species_count_significant <- lapply(combined_sampling_fractions_species_per_genus_significant, function(x) sum(x$Total_Species))
total_species_count_not_significant <- lapply(combined_sampling_fractions_species_per_genus_not_significant, function(x) sum(x$Total_Species))

# Now I can unlist those lists
total_species_count_significant <- unlist(total_species_count_significant)
total_species_count_not_significant <- unlist(total_species_count_not_significant)

# Now I can do a t-test to see if the total number of species in the significant and non-significant orders are different
t.test(total_species_count_significant, total_species_count_not_significant, alternative = "two.sided")

# ggplot boxplot with jitter
# Create a data frame with the total species count and significance
total_species_count <- data.frame(
  Total_Species_Count = c(total_species_count_significant, total_species_count_not_significant),
  Significance = rep(c("Significant", "Non-significant"), c(length(total_species_count_significant), length(total_species_count_not_significant)))
)

# Plot with boxplot and jitter
ggplot(total_species_count, aes(x = Significance, y = Total_Species_Count, fill = Significance)) +
  geom_boxplot(alpha = 0.9, fill = cols) +
  geom_jitter(color = "black", size = 1, alpha = 0.6, width = 0.2) +
  labs(
	   x = "Significance",
	   y = "Total Species Count") +
  theme_minimal() +
  theme(
	legend.position = "none",
	plot.title = element_text(size = 11)
  )

# save the plot
ggsave(file.path(plot_folder, "total_species_count_plot.png"), width = 6, height = 4, dpi = 300)

# Now we do the same test but with species_count
# I can do the same with the species count
species_count_significant <- lapply(combined_sampling_fractions_species_per_genus_significant, function(x) sum(x$Species_Count))
species_count_not_significant <- lapply(combined_sampling_fractions_species_per_genus_not_significant, function(x) sum(x$Species_Count))

# Now I can unlist those lists
species_count_significant <- unlist(species_count_significant)
species_count_not_significant <- unlist(species_count_not_significant)

# Now I can do a t-test to see if the total number of species in the significant and non-significant orders are different
t.test(species_count_significant, species_count_not_significant, alternative = "two.sided")

# ggplot boxplot with jitter
# Create a data frame with the total species count and significance
species_count <- data.frame(
  Species_Count = c(species_count_significant, species_count_not_significant),
  Significance = rep(c("Significant", "Non-significant"), c(length(species_count_significant), length(species_count_not_significant)))
)

# Plot with boxplot and jitter
ggplot(species_count, aes(x = Significance, y = Species_Count, fill = Significance)) +
  geom_boxplot(alpha = 0.9, fill = cols) +
  geom_jitter(color = "black", size = 1, alpha = 0.6, width = 0.2) +
  labs(
	   x = "Significance",
	   y = "Species Count") +
  theme_minimal() +
  theme(
	legend.position = "none",
	plot.title = element_text(size = 11)
  )

# save the plot
ggsave(file.path(plot_folder, "species_count_plot.png"), width = 6, height = 4, dpi = 300)

sort(species_count$Species_Count)
species_count$Species_Count
