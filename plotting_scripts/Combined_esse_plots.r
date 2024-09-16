# Load libraries
library(ggplot2)
library(dplyr)
library(maps)
library(gridExtra)
library(data.table)
library(cowplot)
library(tidyr)
library(ggridges)
library(MetBrewer)
library(ape)

# Get data dir
data_dir <- "/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output/Esse_done_results/"

# Plots directory
plots_dir <- "/home/au543206/GenomeDK/Trf_models/workflow/05_figures/"

# Get a list of all the files in the folder with the pattern combined_esse_runs_", taxon, ".csv"
files_list <- list.files(data_dir, full.names = TRUE, recursive = TRUE, pattern = 'combined_esse_runs_')

# Now we want the taxon names from the file names
taxon <- gsub(".*combined_esse_runs_(.*).csv", "\\1", files_list)
taxon

# Remove the ones which contains _0.33.log
taxon <- taxon[!grepl("_0.33.log", taxon)]
files_list <- files_list[!grepl("_0.33.log", files_list)]

# I now want to load each of these files and combine them into a single data frame
# I want one of the rows in the dataframe to be the taxon name

data <- data.frame()

# Looping through each file in the files_list
for (i in 1:length(files_list)) {
  print(paste("Processing file", i, "of", length(files_list)))
  file <- files_list[i]

  df <- fread(file)
  #print(head(df))
  df$taxon <- taxon[[i]]
  data <- rbind(data, df)
}

# Now I want to create a violin plot of the data
data_lambda <- data %>% dplyr::select(lambda_A_0, lambda_B_0, lambda_W_0, taxon)
data_lambda

data_mu <- data %>% dplyr::select(mu_A_0, mu_B_0, loss_A_0, loss_B_0, taxon)
head(data_mu)

# OKay so i need to melt the data so that I can plot it and what I think I need to do is combine Lambda_A_0, lambda_B_0, Lambda_W_0 into a single column
# Reshape the data to long format
data_lambda_long <- data_lambda %>%
  pivot_longer(cols = starts_with("lambda"), 
               names_to = "lambda", 
               values_to = "value")


# Convert lambda_A_0 to Tropical Rainforest
data_lambda_long$lambda <- gsub("lambda_A_0", "Tropical Rainforest", data_lambda_long$lambda)

# Convert lambda_B_0 to Outside Tropical Rainforest
data_lambda_long$lambda <- gsub("lambda_B_0", "Outside Tropical Rainforest", data_lambda_long$lambda)

# Convert lambda_W_0 to Widespread
data_lambda_long$lambda <- gsub("lambda_W_0", "Widespread", data_lambda_long$lambda)

head(data_lambda_long)

cols <- met.brewer("Java", n=3, type = "discrete")

trfcol <- cols[3]
non_trfcol <- cols[1]
widespread <- cols[2]


# With transparency (right)
# can we make both the x and y axes  log scaled? 

p2 <- ggplot(data=data_lambda_long, aes(x=value, group=lambda, fill=lambda)) +
    geom_density(adjust=1.5, alpha=.7, ) +
	scale_fill_manual(values = c(non_trfcol,trfcol, widespread)) +
    theme_minimal() +
	scale_x_continuous(limits = c(0,25), expand = c(0,0)) +
	theme(  legend.position = "bottom",
			axis.text.y = element_blank(),
			axis.title.y = element_blank()
			) +
	xlab("Speciation Rate") + 
	labs(fill = "")
print(p2)

p2_log <- ggplot(data_lambda_long, aes(x = value, group = lambda, fill = lambda))+
    geom_density(adjust=1.5, alpha=.7, ) +
	scale_fill_manual(values = c(non_trfcol,trfcol, widespread)) +
    theme_minimal() +
	scale_x_log10() +
	theme(legend.position = "bottom")
print(p2_log)


# Yes speciation in Trf and widespread are higher than non_trf and 0
# There are a lot of residual variation suggesting that there are plenty other factors influencing speciation rates
# The variability in in intercepts across the taxons suggests that each taxon differ in their baseline speciation rates and supports the use of them as random effects. 


# Now I want to try and create the same plot but with the different taxons
# What I need to do is create a dataset where we have 1 column which is all the values of lambda_A_0, lambda_B_0, lambda_W_0
# THen a second row which is the taxon name and then a third row which is the lambda name
data_for_violin_plot <- data_lambda_long %>% dplyr::select(value, taxon, lambda)
data_for_violin_plot

# Now we want to create violin plots for each of the taxons across the 3 different lambda categories.
# Assuming `data_for_violin_plot` is already in the long format as described earlier
# p3 <- ggplot(data=data_for_violin_plot, aes(x=log(value), y=lambda, fill=lambda)) +
#   geom_violin(adjust=1.5, alpha=.4) +
#   theme_minimal() +
#   scale_x_continuous(limits = c(-7,7)) +
#   facet_wrap(~taxon, scales = "free_y") +
#   labs(x = "Log(Value)", y = "Lambda Category")
# p3



# Plotting
p3_facet <- ggplot(data=data_for_violin_plot, aes(x=lambda, y=log10(value))) +
  geom_violin(aes(fill=lambda), adjust=1.5, alpha=.7, position=position_dodge(width=0.9)) +
  scale_fill_manual(values = c(non_trfcol,trfcol, widespread)) +
  theme_minimal() +
  facet_wrap(~taxon) +
  labs(x = "Lambda Category", y = "Value", title = "Violin Plots of Lambda Values by Taxon") +
  theme(
	legend.position = "none",
	axis.text.x = element_blank()
	) # Remove legend if not needed
p3_facet

#
p3_facet_test <- ggplot(data=data_for_violin_plot, aes(x=lambda, y=log10(value))) +
  geom_violin(aes(fill=lambda), adjust=1.5, alpha=.7, position=position_dodge(width=0.9)) +
  scale_fill_manual(values = c(non_trfcol,trfcol, widespread)) +
  theme_minimal() +
  facet_wrap(~taxon) +
  labs(x = "Lambda Category", y = "Value", title = "Violin Plots of Speciation by Taxon") +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    strip.text = element_text(size = 5) # Adjust facet title size
  )
p3_facet_test

# Save the plot in A4 landscape format
ggsave(paste0(plots_dir,"p3_facet_test_a4_landscape.png"), 
       plot = p3_facet_test, 
       width = 11.69,  # A4 width in inches
       height = 8.27,  # A4 height in inches
       units = "in")

# Now I do the same thing with extinction rates




# p3_stacked <- ggplot(data = data_for_violin_plot,aes(lambda, log(value), identy = taxon, fill = lambda)) + 
#   geom_violin(position = "identity") +
#   theme_minimal() + 
#   theme(legend.position = "none") +
#   coord_flip()
# p3_stacked

# Create a dataframe where we calculate the average value of lambda_A_0, lambda_B_0, lambda_W_0 for each taxon
data_avg_lambda <- data_lambda %>%
  group_by(taxon) %>%
  summarise(avg_lambda_A_0 = mean(lambda_A_0, na.rm = TRUE),
			avg_lambda_B_0 = mean(lambda_B_0, na.rm = TRUE),
			avg_lambda_W_0 = mean(lambda_W_0, na.rm = TRUE))

# Convert to long format
data_avg_lambda_long <- data_avg_lambda %>%
  pivot_longer(cols = starts_with("avg_lambda"), 
			   names_to = "lambda", 
			   values_to = "value")


# p_jitter <- ggplot(data=data_avg_lambda_long, aes(x=lambda, y=value, color=lambda)) +
#   geom_jitter(width = 0.2, size = 3) +
#   theme_minimal() +
#   #facet_wrap(~taxon, scales = "free_y") +
#   labs(x = "Lambda Category", y = "Log(Value)", title = "Average Lambda Values by Taxon") +
#   theme(legend.position = "none") + # Remove legend if not needed
#   coord_flip()
# p_jitter

# Maybe We should try with Ridgeline plots, I can create Ridgeline plots for each of the 3 lambda values and then cowplot them on top of each other.
p_ridge_a <- ggplot(data = data_lambda, aes(x = lambda_A_0, y = taxon)) + 
			geom_density_ridges(scale = 10, fill = trfcol, color = "white") +
			theme(legend.position = "none") +
			theme_minimal() +
			theme(
				axis.text.y = element_blank(),
				axis.title.y = element_blank(),
				axis.text.x = element_blank(),
			) +
			xlab("") +
			scale_y_discrete(expand = expansion(add = c(0,25))) +
			scale_x_continuous(limits = c(0,25), expand = c(0,0))
p_ridge_a

p_ridge_b <- ggplot(data = data_lambda, aes(x = lambda_B_0, y = taxon)) + 
			geom_density_ridges(scale = 10, fill = non_trfcol, colour = "white") +
			theme(legend.position = "none") +
			theme_minimal() +
			theme(
				axis.text.y = element_blank(),
				axis.title.y = element_blank(),
				axis.text.x = element_blank(),
			) +
			xlab("") +
			scale_x_continuous(limits = c(0,25), expand = c(0,0))
p_ridge_b

p_ridge_w <- ggplot(data = data_lambda, aes(x = lambda_W_0, y = taxon)) + 
			geom_density_ridges(scale = 10, fill = widespread, colour = "white") +
			theme(legend.position = "none") +
			theme_minimal() +
			theme(
				axis.text.y = element_blank(),
				axis.title.y = element_blank()
			) + 
			xlab("") +
			scale_x_continuous(limits = c(0,25), expand = c(0,0))
p_ridge_w

# Can I now create a cowplot of these 3 plots where they are all above each other
p_ridge <- plot_grid(p_ridge_a, p_ridge_b, p_ridge_w, ncol = 1)
p_ridge

# Now I want to add p2 to the bottom of p_ridge
p_ridge_dens <- plot_grid(
						p_ridge + theme(plot.background = element_rect(color = "black")),
						p2 + theme(plot.background = element_rect(color = "black")),
						ncol = 1, labels = c("A", "B"))
p_ridge_dens

# Save the plot as a pdf file
#pdf(paste0(plots_dir,"combined_esse_plots.pdf"), width = 5, height = 10)
png(paste0(plots_dir,"combined_esse_plots.png"), width = 600, height = 1200)
p_ridge_dens

dev.off()


# # Maybe We should try with Ridgeline plots, I can create Ridgeline plots for each of the 3 lambda values and then cowplot them on top of each other.
# p_ridge_a_log <- ggplot(data = data_lambda, aes(x = lambda_A_0, y = taxon)) + 
# 			geom_density_ridges(scale = 10, fill = trfcol, color = "white") +
# 			theme(legend.position = "none") +
# 			theme_minimal() +
# 			theme(
# 				axis.text.y = element_blank(),
# 				axis.title.y = element_blank(),
# 				axis.text.x = element_blank(),
# 			) +
# 			xlab("") +
# 			scale_x_log10( expand = c(0,0))
# p_ridge_a_log


# p_ridge_b_log <- ggplot(data = data_lambda, aes(x = lambda_B_0, y = taxon)) + 
# 			geom_density_ridges(scale = 10, fill = non_trfcol, colour = "white") +
# 			theme(legend.position = "none") +
# 			theme_minimal() +
# 			theme(
# 				axis.text.y = element_blank(),
# 				axis.title.y = element_blank(),
# 				axis.text.x = element_blank(),
# 			) +
# 			xlab("") +
# 			scale_x_log10(expand = c(0,0))
# p_ridge_b_log

# p_ridge_w_log <- ggplot(data = data_lambda, aes(x = lambda_W_0, y = taxon)) + 
# 			geom_density_ridges(scale = 10, fill = widespread, colour = "white") +
# 			theme(legend.position = "none") +
# 			theme_minimal() +
# 			theme(
# 				axis.text.y = element_blank(),
# 				axis.title.y = element_blank(),
# 				axis.text.x = element_blank(),
# 			) + 
# 			scale_x_log10(expand = c(0,0))
# p_ridge_w_log

# # Can I now create a cowplot of these 3 plots where they are all above each other
# p_ridge_log <- plot_grid(p_ridge_a_log, p_ridge_b_log, p_ridge_w_log, ncol = 1)
# p_ridge_log


# # Now make a 
# p_ridge_dens_log <- plot_grid(
# 								p_ridge_log + theme(plot.background = element_rect(color = "black")),
# 								p2_log + theme(plot.background = element_rect(color = "black")),
# 								ncol = 1, labels = c("A", "B"))
# p_ridge_dens_log

# # Now I want to add p2 to the bottom of p_ridge
# p_ridge_dens <- plot_grid(
# 						p_ridge + theme(plot.background = element_rect(color = "black")),
# 						p2 + theme(plot.background = element_rect(color = "black")),
# 						ncol = 1, labels = c("A", "B"))
# p_ridge_dens


# p_ridge_dens_log




# # Comparing plots
# p_ridge_dens
# p_ridge_dens_log


#########################################################################################################################
################################################## Plotting Mu data ####################################################
#########################################################################################################################

# Reshape the data to long format
data_mu_long <- data_mu %>%
  pivot_longer(cols = starts_with("mu"), 
			   names_to = "mu", 
			   values_to = "value")

# Convert mu_A_0 to Tropical Rainforest
data_mu_long$mu <- gsub("mu_A_0", "Tropical Rainforest", data_mu_long$mu)

# Convert mu_B_0 to Outside Tropical Rainforest
data_mu_long$mu <- gsub("mu_B_0", "Outside Tropical Rainforest", data_mu_long$mu)

#  Now make a density plot
p4 <- ggplot(data=data_mu_long, aes(x=value, group=mu, fill=mu)) +
	geom_density(adjust=1.5, alpha=.7, ) +
	scale_fill_manual(values = c(non_trfcol,trfcol)) +
	theme_minimal() +
	scale_x_continuous(limits = c(0,25), expand = c(0,0)) +
	theme(  legend.position = "bottom",
			axis.text.y = element_blank(),
			axis.title.y = element_blank()
			) +
	xlab("Extinction Rate") + 
	labs(fill = "")

p4

# Now I want to create the ridgeline plots for the extinction rates
p_ridge_mu_a <- ggplot(data = data_mu, aes(x = mu_A_0, y = taxon)) + 
			geom_density_ridges(scale = 10, fill = trfcol, color = "white") +
			theme(legend.position = "none") +
			theme_minimal() +
			theme(
				axis.text.y = element_blank(),
				axis.title.y = element_blank(),
			) +
			xlab("") +
			scale_y_discrete(expand = expansion(add = c(0,25))) +
			scale_x_continuous(limits = c(0,25), expand = c(0,0))
p_ridge_mu_a

p_ridge_mu_b <- ggplot(data = data_mu, aes(x = mu_B_0, y = taxon)) + 
			geom_density_ridges(scale = 10, fill = non_trfcol, colour = "white") +
			theme(legend.position = "none") +
			theme_minimal() +
			theme(
				axis.text.y = element_blank(),
				axis.title.y = element_blank(),
			) +
			xlab("") +
			scale_x_continuous(limits = c(0,25), expand = c(0,0))

p_ridge_mu_b


# Can I now create a cowplot of these 3 plots where they are all above each other
p_ridge_mu <- plot_grid(p_ridge_mu_a, p_ridge_mu_b, ncol = 1)

# Adding the density plot to the bottom of the ridgeline plots
p_ridge_mu_dens <- plot_grid(
						p_ridge_mu + theme(plot.background = element_rect(color = "black")),
						p4 + theme(plot.background = element_rect(color = "black")),
						ncol = 1, labels = c("C", "D"))

p_ridge_mu_dens

# save the plot as a png file
png(paste0(plots_dir,"combined_mu_plots.png"), width = 600, height = 1200)
p_ridge_mu_dens
dev.off()

# Now combine the plots of speciation and extinction in 2 columns
p_combined <- plot_grid(p_ridge_dens, p_ridge_mu_dens, ncol = 2)
p_combined

# Save the plot as a png file in A4 format
png(paste0(plots_dir,"combined_esse_mu_plots.png"), width = 900, height = 1200)
p_combined
dev.off()


###################################################################################################
###################################################################################################
###################################################################################################
# Testing wether there is a significant effect of 
# Can we test if there is a significant difference between the lambda values for the different taxons
library(lmerTest)

head(data_lambda_long)

model_lmer <- lmer(value ~ lambda + (1|taxon), data = data_lambda_long)
summary(model_lmer) # Does this model actually suggest that Trf and widespread have higher speciation rates than non_trf?

# I need to select a random species from each subphylogeny i.e the unique names in the taxon column
# Then I need to prune the Smith and brown phylogeny to only include these species.
# I will then rename each of the tips to fit the taxon label in data_lambda_long

# First I need to find a random species from each subphylogeny
# I think the easiest way to do this is to find the unique "taxons" in the data_lambda_long dataframe
# Then find the Esse trees and which contains the taxon in its name.
# Then I select a random species from that tree and add it to a list of species to keep along with its taxon 
# I then prune the Smith and Brown phylogeny to only include these species and then rename the tips to the taxon name

# Find list of all the unique taxons
taxons <- unique(data_lambda_long$taxon)

# Load the Esse data
# Esse on orders
Esse_order_folder <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/orders"

# Esse on families
Esse_family_folder <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/families"

# Esse on subfamilies
Esse_subfamily_folder <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/subset_of_orders"

# The pattern we want to find in each folder is 
pattern <- "_Esse_tree.tre"

# Get a list of all Rdata files in the folder
file_list_1 <- list.files(Esse_order_folder, pattern = pattern, full.names = TRUE)
file_list_2 <- list.files(Esse_family_folder, pattern = pattern, full.names = TRUE)
file_list_3 <- list.files(Esse_subfamily_folder, pattern = pattern, full.names = TRUE)

file_list_esse_trees <- c(file_list_1, file_list_2, file_list_3)
file_list_esse_trees
# Find all trees in file list which contains the names in taxons
trees_which_i_need <- list()
for (taxon in taxons) {
  tree <- grep(paste0(taxon,"_Esse_tree.tre"), file_list_esse_trees, value = TRUE)
  trees_which_i_need <- c(trees_which_i_need, tree)
}
trees_which_i_need

# Find the species for each tree and save the species and the taxon in a list
species_to_keep <- c()
taxon_list <- c()
for (tree in trees_which_i_need) {
  tree_read <- read.tree(tree)
  species <- sample(tree_read$tip.label, 1)
  print(species)
  
  split_filepath <- strsplit(tree, "/")
  last_element <- tail(split_filepath[[1]], n = 1)
  taxon <- gsub("_Esse_tree.tre", "", last_element)

  species_to_keep <- c(species_to_keep, species)
  taxon_list <- c(taxon_list, taxon)
}
tip_taxon_df <- data.frame(species = species_to_keep, taxon = taxon_list)
tip_taxon_df

# Loop through all the unique taxons in tip_taxon_df and only keep one of the rows with each taxon
for (taxon in unique(tip_taxon_df$taxon)) {
  
  # Find the rows which contain the taxon
  rows <- which(tip_taxon_df$taxon == taxon)
  print(paste0(taxon, " ", length(rows)))

  # Select a random row
  if (length(rows) > 1) {

  random_row <- sample(rows, 1)
  rows <- rows[-which(rows == random_row)]
  print(random_row)
  print(rows)
  # Remove all the other rows
  tip_taxon_df <- tip_taxon_df[-rows,]
  } else {
	next
  }
}
tip_taxon_df

# Load the renamed Smith and Brown phylogeny
smb_tree <- read.tree("/home/au543206/GenomeDK/Trf_models/data/GBMB_pruned.tre")

# Now I want to prune the Smith and Brown phylogeny to only include the species in species_to_keep
smb_tree_pruned <- drop.tip(smb_tree, setdiff(smb_tree$tip.label, tip_taxon_df$species))

# Match the species names in the tree with those in the dataframe
match_idx <- match(smb_tree_pruned$tip.label, tip_taxon_df$species)

# Replace the tip labels in the tree with the corresponding taxon names
smb_tree_pruned$tip.label <- tip_taxon_df$taxon[match_idx]

# Calculate covariance matrix for the random variable.
phylo_cov_matrix <- vcv.phylo(smb_tree_pruned, corr = TRUE)

rownames(phylo_cov_matrix)

# Do all the unique taxa have a corresponding value in the phylo_cov_matrix
all(data_avg_lambda_long$taxon %in% rownames(phylo_cov_matrix))

# Which data_avg_lambda_long$taxon are not in the rownames of phylo_cov_matrix
missing_from_phyl <- data_avg_lambda_long$taxon[!data_avg_lambda_long$taxon %in% rownames(phylo_cov_matrix)]
missing_from_phyl

# Fucking boot them from the dataframe I dont care ;(
data_avg_lambda_long <- data_avg_lambda_long[-which(data_avg_lambda_long$taxon %in% missing_from_phyl),]

# Do all the unique taxa have a corresponding value in the phylo_cov_matrix
all(data_avg_lambda_long$taxon %in% rownames(phylo_cov_matrix))

# Do i need to select the mean value of lambda for each taxon as the MCMCglmm might not like that there are multiple values for each taxon
# I think I need to calculate the mean value of lambda for each taxon
data_lambda_long_mean <- data_lambda_long %>%
  group_by(taxon, lambda) %>%
  summarise(mean_value = mean(value, na.rm = TRUE))

# Setting taxon to be a factor
data_lambda_long_mean$taxon <- as.factor(data_lambda_long_mean$taxon)

phylo_cov_matrix <- phylo_cov_matrix[data_lambda_long_mean$taxon, data_lambda_long_mean$taxon]

# Load phyr
library(phyr)

# Fit a phylogenetic mixed model with random effects
model_phyr <- pglmm(mean_value ~ lambda + (1|taxon), 
                    data = data_lambda_long_mean, 
                    cov_ranef = list(taxon = smb_tree_pruned),  # Provide the phylogenetic tree
                    family = "gaussian")

# Check summary of the model
summary(model_phyr)


# Now I want to test if there is a significant difference between the mu values for the different taxons
data_mu_long_mean <- data_mu_long %>%
  group_by(taxon, mu) %>%
  summarise(mean_value = mean(value, na.rm = TRUE))


# Setting taxon to be a factor
data_mu_long_mean$taxon <- as.factor(data_mu_long_mean$taxon)

phylo_cov_matrix <- phylo_cov_matrix[data_mu_long_mean$taxon, data_mu_long_mean$taxon]

# Fit a phylogenetic mixed model with random effects
model_phyr_mu <- pglmm(mean_value ~ mu + (1|taxon), 
					data = data_mu_long_mean, 
					cov_ranef = list(taxon = smb_tree_pruned),  # Provide the phylogenetic tree
					family = "gaussian")

summary(model_phyr_mu)

# Convert the p values to non scientific notation
model_phyr[9] <- format(model_phyr[9], scientific=FALSE)

model_phyr[9][[1]][1]
model_phyr[9][[1]][2]
model_phyr[9][[1]][3]

format(max_speciation_rate, scientific=FALSE) #343.5 

format(model_phyr[9][[1]][1], scientific=FALSE)
format(model_phyr[9][[1]][2], scientific=FALSE)
format(model_phyr[9][[1]][3], scientific=FALSE)
