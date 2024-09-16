# Load necessary libraries
library(ggplot2)

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

# Can we for each of these taxon names run a model which tests to see if there is a significant difference between the lambda and mu values?
# I think we can do this by using a t-test
# We can use the t.test function in R to do this
# The t.test function requires a numeric vector as input

# For each of the taxon names we want to run 2 t-tests one for lambda and one for mu
# We want to test whether there is a difference between lambda_A_0, lambda_B_0 and lambda_W_0.
# We also want to test whether there is a difference between mu_A_0 and mu_B_0.
list_of_aov_results <- list()

for(i in 1:length(unique(data$taxon))) {
  #print progress
  print(paste("Processing taxon", i, "of", length(unique(data$taxon))))

  # Get the data for the taxon
  taxon_data <- data[data$taxon == unique(data$taxon)[i],]
  
  # Get the lambda values
  lambda_A_0 <- taxon_data$lambda_A_0
  lambda_B_0 <- taxon_data$lambda_B_0
  lambda_W_0 <- taxon_data$lambda_W_0
  
  # Get the mu values
  mu_A_0 <- taxon_data$mu_A_0
  mu_B_0 <- taxon_data$mu_B_0
  
  # Test whether the data is normally distributed
  lambda_A_0_norm <- ks.test(lambda_A_0, "pnorm", mean(lambda_A_0), sd(lambda_A_0))	
  lambda_B_0_norm <- ks.test(lambda_B_0, "pnorm", mean(lambda_B_0), sd(lambda_B_0))
  lambda_W_0_norm <- ks.test(lambda_W_0, "pnorm", mean(lambda_W_0), sd(lambda_W_0))

  #making long lambda data
  long_data <- melt(taxon_data[, .(lambda_A_0, lambda_B_0, lambda_W_0)], variable.name = "group", value.name = "lambda")

  # Test to see if they are all normally distributed
  if(all(lambda_A_0_norm$p.value > 0.05, lambda_B_0_norm$p.value > 0.05, lambda_W_0_norm$p.value > 0.05)){ # Normal data

	# Run the levene test of equal variance
	if(leveneTest(lambda ~ group, data = long_data)$p.value > 0.05){ # Normal data and equal variance
	# Normal anova
	anova_result_lambda <- aov(lambda ~ group, data = long_data)

  } else { # Normal data but not equal variance

	# Welch's anmova
	#anova_result_lambda <- oneway.test(lambda ~ group, data = long_data, var.equal = FALSE)
	anova_result_lambda <- aov(lambda ~ group, data = long_data, var.equal = FALSE)
  		}

	} else { # Non-normal data 
	#anova_result_lambda <- oneway.test(lambda ~ group, data = long_data, var.equal = FALSE)
	anova_result_lambda <- aov(lambda ~ group, data = long_data, var.equal = FALSE)
	}

# Do the same for mu values
mu_A_0_norm <- ks.test(mu_A_0, "pnorm", mean(mu_A_0), sd(mu_A_0))
mu_B_0_norm <- ks.test(mu_B_0, "pnorm", mean(mu_B_0), sd(mu_B_0))

#making long mu data
long_data <- melt(taxon_data[, .(mu_A_0, mu_B_0)], variable.name = "group", value.name = "mu")

# Test to see if they are all normally distributed
if(all(mu_A_0_norm$p.value > 0.05, mu_B_0_norm$p.value > 0.05)){ # Normal data

	# Run the levene test of equal variance
	if(leveneTest(mu ~ group, data = long_data)$p.value > 0.05){ # Normal data and equal variance
	# Normal anova
	anova_result_mu <- aov(mu ~ group, data = long_data)

  } else { # Normal data but not equal variance

	# Welch's anmova
	#anova_result_mu <- oneway.test(mu ~ group, data = long_data, var.equal = FALSE)
	anova_result_mu <- aov(mu ~ group, data = long_data, var.equal = FALSE)
  		}

	} else { # Non-normal data 
	#anova_result_mu <- oneway.test(mu ~ group, data = long_data, var.equal = FALSE)
	anova_result_mu <- aov(mu ~ group, data = long_data, var.equal = FALSE)
	}

# Save the results to a list
list_of_aov_results[[unique(data$taxon)[i]]] <- list(anova_result_lambda, anova_result_mu)
}

list_of_aov_results

# How many phylogenies have significant differences between lambda_A_0, lambda_B_0 and lambda_W_0?
# How many phylogenies have significant differences between mu_A_0 and mu_B_0?

# Initialize a counter
lambda_significant <- 0
mu_significant <- 0

# Loop through the list of results
for (i in 1:length(list_of_aov_results)) {
  # Get the results for the lambda values
  lambda_results <- list_of_aov_results[[i]][[1]]
  
  # Get the results for the mu values
  mu_results <- list_of_aov_results[[i]][[2]]
  
  # Check if the lambda results are significant
  if (lambda_results$p.value < 0.05) {
	lambda_significant <- lambda_significant + 1
  }
  
  # Check if the mu results are significant
  if (mu_results$p.value < 0.05) {
	mu_significant <- mu_significant + 1
  }
}

lambda_significant
mu_significant

#


# find the max speciation rate
max_speciation_rate <- max(data$lambda_A_0, data$lambda_B_0, data$lambda_W_0)
format(max_speciation_rate, scientific=FALSE) #343.5 

#find the min speciation rate
min_speciation_rate <- min(data$lambda_A_0, data$lambda_B_0, data$lambda_W_0)
format(min_speciation_rate, scientific=FALSE) #0.000000

#find the max extinction rate
max_extinction_rate <- max(data$mu_A_0, data$mu_B_0)
format(max_extinction_rate, scientific=FALSE) #0.000000

#find the min extinction rate
min_extinction_rate <- min(data$mu_A_0, data$mu_B_0)
format(min_extinction_rate, scientific=FALSE) #0.000000

# Now we need to do post hoc tests to see which of the groups are significantly different from each other
# We can do this by using the TukeyHSD function in R

# Create a list to store the results
list_of_tukey_results <- list()

# Loop through the list of results
for (i in 1:length(list_of_aov_results)) {
  # Print progress
  print(paste("Processing taxon", i, "of", length(list_of_aov_results)))

  # Get the results for the lambda values
  lambda_results <- list_of_aov_results[[i]][[1]]
  
  # Get the results for the mu values
  mu_results <- list_of_aov_results[[i]][[2]]
  
  # Check if the lambda results are significant
  if (summary(lambda_results)[[1]][[5]][[1]] < 0.05) {
	# Run the TukeyHSD test
	tukey_results_lambda <- TukeyHSD(lambda_results)
	
	# Save the results to a list
	#list_of_tukey_results[[unique(data$taxon)[i]]] <- tukey_results_lambda
  }
  
  # Check if the mu results are significant
  if (summary(mu_results)[[1]][[5]][[1]] < 0.05) {
	# Run the TukeyHSD test
	tukey_results_mu <- TukeyHSD(mu_results)
  }

  if(summary(mu_results)[[1]][[5]][[1]] < 0.05 & summary(lambda_results)[[1]][[5]][[1]] < 0.05){
  list_of_tukey_results[[unique(data$taxon)[i]]] <- list(tukey_results_lambda, tukey_results_mu)
  } else if (summary(mu_results)[[1]][[5]][[1]] < 0.05 & summary(lambda_results)[[1]][[5]][[1]] > 0.05){
	 list_of_tukey_results[[unique(data$taxon)[i]]] <- tukey_results_mu
  } else if (summary(mu_results)[[1]][[5]][[1]] > 0.05 & summary(lambda_results)[[1]][[5]][[1]] < 0.05){
	 list_of_tukey_results[[unique(data$taxon)[i]]] <- tukey_results_lambda
  }
}


Trf_higher_lambda <- 0
Trf_higher_mu <- 0

non_trf_higher_lambda <- 0
non_trf_higher_mu <- 0

# Loop through the list of results and count the number of times Trf has higher lambda and mu values
for (i in 1:length(list_of_tukey_results)) {
  # Get the results for the taxon
  results <- list_of_tukey_results[[i]]
  
  # Check if the results are a list
  if (is.list(results)) {
	# Get the results for the lambda values
	tukey_results_lambda <- results[[1]]
	
	# Get the results for the mu values
	tukey_results_mu <- results[[2]]
	
	# Check if the Trf group has higher lambda values
	if (tukey_results_lambda$group[1] < 0) {
	  Trf_higher_lambda <- Trf_higher_lambda + 1
	} else {
	  non_trf_higher_lambda <- non_trf_higher_lambda + 1
	}
	
	# Check if the Trf group has higher mu values
	if (tukey_results_mu$group[1] < 0) {
	  Trf_higher_mu <- Trf_higher_mu + 1
	} else {
	  non_trf_higher_mu <- non_trf_higher_mu + 1
	}
  } else {
	
	# Get the i value of the results
	print(paste0("i: ", i))

	# Find out of the element is of lambda or mu
	if (grepl("lambda", names(results))) {
	  # Get the results for the lambda values
	  tukey_results_lambda <- results
	  
	  # Check if the Trf group has higher lambda values
	  if (tukey_results_lambda$group[1] < 0) {
		Trf_higher_lambda <- Trf_higher_lambda + 1
	  } else {
		non_trf_higher_lambda <- non_trf_higher_lambda + 1
	  }
	} else {
	  # Get the results for the mu values
	  tukey_results_mu <- results
	  
	  # Check if the Trf group has higher mu values
	  if (tukey_results_mu$group[1] < 0) {
		Trf_higher_mu <- Trf_higher_mu + 1
	  } else {
		non_trf_higher_mu <- non_trf_higher_mu + 1
	  }
	}
  }

	# Get the results for the mu values
	tukey_results_mu <- results

  }
}


Trf_higher_lambda
Trf_higher_mu 

non_trf_higher_lambda
non_trf_higher_mu

