#########################################################################################################################
################################################## Loading packages  ####################################################
#########################################################################################################################

# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("data.table", "coda")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

#########################################################################################################################
################################################## Loading the data #####################################################
#########################################################################################################################

# specify the file paths of the Esse runs
data_dir <- "/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output" # When running it on local machine

# Get a list of all the files in the folder
# We need to use regex to find all the files that end in a number and then .log
esse_runs <- list.files(data_dir, pattern = "*_[0-9]_hidden_states_0.33_[0-9].log", full.names = TRUE)
esse_runs

# Now we need to find all the familes/subfamilies/orders
# This is just the part that is between Esse_states_ and _0.33_[0-9].log
# We can use regex to extract this information
esse_runs_info <- unique(gsub(".*Esse_output_|_0.33_[0-9].log", "", esse_runs))
esse_runs_info

# Can we remove the hidden states from the esse_runs_info
esse_runs_names <- gsub("_hidden_states", "", esse_runs_info)
esse_runs_names

# For each element in esse_runs_info, we need to find the matching files in esse_runs
# We can use grep to do this
esse_runs_list <- lapply(esse_runs_info, function(x) {
  grep(x, esse_runs, value = TRUE)
})

# For now we need to remove the ones without hidden_states in their name.
esse_runs_list <- lapply(esse_runs_list, function(x) {
  x[grep("hidden_states", x)]
})

# for now lets remove elements from esse_run_list which have a lenght shorter than 2.
esse_runs_list <- lapply(esse_runs_list, function(x) {
  if (length(x) > 1) {
	return(x)
  }
})

esse_runs_list

# Now for each of the files found in each list element of esse_runs_list, we need to read the data.
# Then we check whether the runs have converged or not.
# We probably need to check all the different combinations of the runs to find the best combination which provides us with the most amount of runs that have converged.
# We can use the gelman.diag function from the coda package to check for convergence.

data <- list()
# Looping through each order/family/subfamily
for (i in 1:length(esse_runs_list)) {
	print(paste(i,"/",length(esse_runs_list)))
	# Looping through each of the runs for the current order/family/subfamily
	esse_fam_list <- list()

	for (j in 1:length(esse_runs_list[[i]])) {
		
		# Checking if is null
		if (is.null(esse_runs_list[[i]][[j]])) {
			next
		}

  		# Read the data
  		data_raw <- read.table(esse_runs_list[[i]][[j]], header = TRUE)

		# data as mcmc object
		data_mcmc <- as.mcmc(data_raw)

		# Append the data to the list
		esse_fam_list[[j]] <- data_mcmc[]

	}
	# Convert to mcmc list
	esse_fam_list <- mcmc.list(esse_fam_list)
	data[[i]] <- esse_fam_list
}


# Can we now calculate the Gelman score for each of the runs?
# If the Gelman score is less than 1.1, we can combine the runs into a single file.
# If the Gelman score is greater than 1.1, we can display a warning message.
# We can use the gelman.diag function from the coda package to calculate the Gelman score.
did_converge <- c()
did_not_converge <- c()
for (i in 1:length(data)){
	print(paste(i,"/",length(data)))

	if (length(data[[i]]) < 2) {
		next
	}
	
	# Calculating the gelman score
	gelman_score <- gelman.diag(data[[i]], multivariate = FALSE)
	
	# Testing to see if they are below 1.1
	if (is.na(gelman_score$psrf[[12,2]]) | is.na(gelman_score$psrf[[13,2]])) {
		print(paste("Esse runs for", esse_runs_info[i], "have no value for beta_lambda."))
		next
	}
	if (gelman_score$psrf[[12,2]] < 1.1 & gelman_score$psrf[[13,2]] < 1.1) {
		combined_data <- data.frame()
		for (j in 1:length(data[[i]])) {
			combined_data <- rbind(combined_data, data[[i]][[j]])
		}
		write.csv(combined_data, paste0("/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output/Esse_done_results/combined_esse_runs_", esse_runs_names[[i]], ".csv"), row.names = FALSE)
		did_converge <- c(did_converge, esse_runs_info[i])

		

	} else {
		# If not converged, display a warning message
		#print(paste("Esse runs for", esse_runs_info[i], "have not converged."))
		did_not_converge <- c(did_not_converge, esse_runs_info[i])
		print(esse_runs_names[[i]])
		print(paste(gelman_score$psrf[[12,2]],gelman_score$psrf[[13,2]]))
	}
}

did_converge
did_not_converge



