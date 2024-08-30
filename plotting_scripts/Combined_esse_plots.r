# Load libraries
library(ggplot2)
library(dplyr)
library(maps)
library(gridExtra)
library(data.table)
library(cowplot)
library(tidyr)
library(ggridges)

# Get data dir
data_dir <- "/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output/Esse_done_results/"

# Plots directory
plots_dir <- "/home/au543206/GenomeDK/Trf_models/workflow/05_figures/"

# Get a list of all the files in the folder with the pattern combined_esse_runs_", taxon, ".csv"
files_list <- list.files(data_dir, full.names = TRUE, recursive = TRUE, pattern = 'combined_esse_runs_')

# Now we want the taxon names from the file names
taxon <- gsub(".*combined_esse_runs_(.*).csv", "\\1", files_list)
taxon

# I now want to load each of these files and combine them into a single data frame
# I want one of the rows in the dataframe to be the taxon name

data <- data.frame()

# Looping through each file in the files_list
for (i in 1:length(files_list)) {
  file <- files_list[i]

  df <- fread(file)
  #print(head(df))
  df$taxon <- taxon[[i]]
  data <- rbind(data, df)
}

# Now I want to create a violin plot of the data from Rubiaceae_3
data_lambda <- data %>% select(lambda_A_0, lambda_B_0, lambda_W_0, taxon)
data_lambda

data_mu <- data %>% select(mu_A_0, mu_B_0, loss_A_0, loss_B_0, taxon)
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

# Can we test if there is a significant difference between the lambda values for the different taxons
library(lmerTest)
model <- lmer(value ~ lambda + (1|taxon), data = data_lambda_long)
summary(model) # Does this model actually suggest that Trf and widespread have higher speciation rates than non_trf?
# Yes speciation in Trf and widespread are higher than non_trf and 0
# There are a lot of residual variation suggesting that there are plenty other factors influencing speciation rates
# The variability in in intercepts across the taxons suggests that each taxon differ in their baseline speciation rates and supports the use of them as random effects. 


# Now I want to try and create the same plot but with the different taxons
# What I need to do is create a dataset where we have 1 column which is all the values of lambda_A_0, lambda_B_0, lambda_W_0
# THen a second row which is the taxon name and then a third row which is the lambda name
data_for_violin_plot <- data_lambda_long %>% select(value, taxon, lambda)
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
	scale_x_continuous(limits = c(0,0.5), expand = c(0,0)) +
	theme(  legend.position = "bottom",
			axis.text.y = element_blank(),
			axis.title.y = element_blank()
			) +
	xlab("Extinction Rate") + 
	labs(fill = "")

p4

# Can we test if there is a significant difference between the mu values for the different taxons
model_mu <- lmer(value ~ mu + (1|taxon), data = data_mu_long)
summary(model_mu) # Does this model actually suggest that Trf and widespread have higher speciation rates than non_trf?
data_mu_long
