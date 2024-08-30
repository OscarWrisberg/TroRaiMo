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

# Establish a list to store the data
data <- list()

# Looping through each order/family/subfamily
for (i in 1:length(esse_runs_list)) {
	print(paste(i, "/", length(esse_runs_list)))

	# Looping through each of the runs for the current order/family/subfamily
	esse_fam_list <- list()

	for (j in 1:length(esse_runs_list[[i]])) {

		# Checking if it is null
		if (is.null(esse_runs_list[[i]][[j]])) {
			next
		}

		# Checking the length of esse_runs_list[[i]]
		if (length(esse_runs_list[[i]]) != 5 && length(esse_runs_list[[i]]) != 10) {
			print("Not the right number of runs")
			next
		}

		# Read the data
		data_raw <- read.table(esse_runs_list[[i]][[j]], header = TRUE)

		# Convert data to mcmc object
		data_mcmc <- as.mcmc(data_raw)

		# Append the data to the list
		esse_fam_list[[j]] <- data_mcmc
	}

	# Check if esse_fam_list has consistent mcmc objects
	if (length(esse_fam_list) > 0) {
		tryCatch({
			esse_fam_list <- mcmc.list(esse_fam_list)
			data[[i]] <- esse_fam_list
		}, error = function(e) {
			print(paste("Error combining MCMC chains for index", i, ": ", e$message))
		})
	}
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
		print(paste("Esse runs for", esse_runs_info[i], "have less than 2 runs."))
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
		did_converge <- c(did_converge, esse_runs_names[i])

	} else {
		# If not converged, display a warning message
		#print(paste("Esse runs for", esse_runs_info[i], "have not converged."))
		did_not_converge <- c(did_not_converge, esse_runs_names[i])
		print(esse_runs_names[[i]])
		print(paste(gelman_score$psrf[[12,2]],gelman_score$psrf[[13,2]]))

	}
}

# Print the number of runs that did and did not converge
did_converge # 88 on 30/08/2024
did_not_converge # 45 on 30/08/2024


####################################################################################

# For the runs that did converge, can we combine the runs from each file into a single file?
# We can use the rbind function to do this.
# We can then write the data to a new file.
for(i in 1:length(did_converge)){
	#print(paste(i,"/",length(did_converge), ":", did_converge[i]))

	# Finding all the files for the current order which converged
	esse_runs_order <- grep(did_converge[i], esse_runs, value = TRUE)

	# Print the esse_runs_order
	print(esse_runs_order)

	# Establish a data frame to store the combined data
	combined_data <- data.frame()

	# Reading the data
	for(j in 1:length(esse_runs_order)){
		data_raw <- read.table(esse_runs_order[j], header = TRUE)
		data_mcmc <- as.mcmc(data_raw)
		print(j)
		if(j == 1){
			combined_data <- data.frame(data_mcmc[])

		} else {
			# Here I need to add some extra to make the column iteration properly go from 0 to x
			# First I need to find the number of iterations in the final row of the combined_data
			combined_data_iteration <- combined_data$Iteration
			max_iteration <- max(combined_data_iteration)

			# Then I will add this number to the iteration column of the data_mcmc
			data_mcmc[,1] <- data_mcmc[,1] + max_iteration

			combined_data <- rbind(combined_data, data_mcmc[])
		}
	}
	# Writing the data to a new file
	name_for_file <- gsub("_hidden_states", "", did_converge[i])
	write.csv(combined_data, paste0("/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output/Esse_done_results/combined_esse_runs_", name_for_file, ".csv"), row.names = FALSE)
}

# Can we now calculate the ESS scores for each of these runs?
# We can use the effectiveSize function from the coda package to calculate the ESS scores.
for (i in 1:length(did_converge)){
	print(paste(i,"/",length(did_converge), ":", did_converge[i]))

	# Reading the data
	data_raw <- read.csv(paste0("/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output/Esse_done_results/combined_esse_runs_", did_converge[i], ".csv"))

	# print the data_raw
	#print(data_raw)

	# Converting the data to a mcmc object
	data_mcmc <- as.mcmc(data_raw)

	# Calculating the ESS scores
	ess_score <- effectiveSize(data_mcmc)

	# Remove iterations
	ess_score <- ess_score[-1]
	
	if(any(ess_score < 200 | is.na(ess_score))){
		print(ess_score[which(ess_score < 200 | is.na(ess_score))])
	}

	# Displaying the ESS scores
	#print(paste0(did_converge,": ",ess_score))
}

#########################################################################################################################
###################################--- Now we deal with the runs that did not converge ---################################
#########################################################################################################################

# Can we make a for loop which checks the convergence for all the combinations of the runs?
# Function to check convergence for a combination of runs
check_convergence <- function(combination, data, threshold) {
  combined_data <- do.call(rbind, lapply(combination, function(index) data[[index]]))
  
  if (any(is.na(combined_data))) {
    print("NA values detected in combined data")
    return(FALSE)
  }
  
  gelman_list <- as.mcmc.list(lapply(combination, function(index) as.mcmc(data[[index]])))
  
  # Check for NA in the mcmc.list
  if (any(sapply(gelman_list, function(x) any(is.na(x))))) {
    print("NA values detected in mcmc.list")
    return(FALSE)
  }
  
  gelman_score <- tryCatch(
    gelman.diag(gelman_list, multivariate = FALSE),
    error = function(e) {
      print(paste("Error in gelman.diag:", e))
      return(NULL)
    }
  )
  
  if (is.null(gelman_score)) {
    return(FALSE)
  }

  # Remove the first row of the psrf matrix
  gelman_psrf <- gelman_score$psrf[-1, ]
  
  if (any(is.na(gelman_psrf))) {
    print(gelman_psrf)
    print("NA values detected in Gelman-Rubin diagnostic")
    return(FALSE)
  }
  
  return(all(gelman_psrf < threshold))
}

# For loop to check the convergence for all the combinations of the runs
list_of_combinations <- list()
for (i in 1:length(did_not_converge)) {
  print(paste(i, "/", length(did_not_converge), ":", did_not_converge[i]))
  
  # Initialize variables to store the best combination and its convergence status
  best_combination <- NULL
  max_converged_runs <- 0
  
  # Load all the files for the current order/family/subfamily
  esse_runs_order <- grep(did_not_converge[i], esse_runs, value = TRUE)
  
  # Load the esse_runs_order as a list of mcmc objects
  data <- lapply(esse_runs_order, function(file) {
    #print(paste("Loading file:", file))
    as.mcmc(read.table(file, header = TRUE))
  })
  
  # Check if data was loaded correctly
  if (any(sapply(data, is.null))) {
    print("Error loading data, skipping to next iteration")
    next
  }
  
  # Loop through all combinations of runs
  for (n in 2:length(data)) { # Start from 2 runs to the total number of runs
    combinations <- combn(1:length(data), n, simplify = FALSE)
    for (comb in combinations) {
      if (check_convergence(comb, data, 1.1)) {
        if (length(comb) > max_converged_runs) {
          max_converged_runs <- length(comb)
          best_combination <- comb
          print(paste("Best combination for", esse_runs_info[i], "is", length(best_combination),"runs"))
        }
      }
    }
  }
# Now I want to save the best combination of runs for each order/family/subfamily
  list_of_combinations[[i]] <- best_combination
}

length(list_of_combinations)

# Can we find the number of runs where the best combination is not null?
good_combinations <- did_not_converge[sapply(list_of_combinations, function(x) !is.null(x))]
length(good_combinations)

# Can we find the files where the best combination of the runs is null?
did_not_converge_null_combinations <- did_not_converge[sapply(list_of_combinations, is.null)]
length(did_not_converge_null_combinations)

# Give the best combinations of runs the names of the orders/families/subfamilies
names(list_of_combinations) <- did_not_converge


# Now run through the list_of_combinations and combine the runs for each order/family/subfamily if the best combination have an ESS above 200 for all variables
combinations_with_good_ess <- list()
combinations_with_bad_ess <- list()
for( i in 1:length(list_of_combinations)){
	# Lets remove hidden states from the name
	taxon <- gsub("_hidden_states", "", names(list_of_combinations)[i])
	print(paste(i,"/",length(list_of_combinations), ":", taxon))

	esse_runs_order <- grep(taxon, esse_runs, value = TRUE)

	#print(esse_runs_order)
	#print(list_of_combinations[[i]])

	# Establish a data frame to store the combined data
	combined_data <- data.frame()

	# Checking if the best combination is NULL
	# if it is we need to go to the next iteration
	if(is.null(list_of_combinations[[i]])){
		next
	}

	# Reading the data
	for (j in 1:length(list_of_combinations[[i]])){
		data_raw <- read.table(esse_runs_order[list_of_combinations[[i]][j]], header = TRUE) 
		data_mcmc <- as.mcmc(data_raw)
		if(j == 1){
			combined_data <- data.frame(data_mcmc[])
		} else {
			# Here I need to add some extra to make the column iteration properly go from 0 to x
			# First I need to find the number of iterations in the final row of the combined_data
			combined_data_iteration <- combined_data$Iteration
			max_iteration <- max(combined_data_iteration)

			# Then I will add this number to the iteration column of the data_mcmc
			data_mcmc[,1] <- data_mcmc[,1] + max_iteration

			combined_data <- rbind(combined_data, data_mcmc[])
		}
	}
	# Writing the data to a new file
	#write.csv(combined_data, paste0("/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output/Esse_done_results/combined_esse_runs_", taxon, ".csv"), row.names = FALSE)

	# Before I write the data to a new file, I need to check the ESS scores
	# Converting the data to a mcmc object
	data_mcmc <- as.mcmc(combined_data)

	# Calculating the ESS scores
	ess_score <- effectiveSize(data_mcmc)

	# Remove iterations
	ess_score <- ess_score[-1]

	if(any(ess_score < 200 | is.na(ess_score))){
		#print(ess_score[which(ess_score < 200 | is.na(ess_score))])
		#print(paste0("ESS scores below 200 for", taxon))
		combinations_with_bad_ess[[taxon]] <- list_of_combinations[[i]]
	} else {
		#print(paste0("All ESS scores above 200 for", taxon))
		combinations_with_good_ess[[taxon]] <- list_of_combinations[[i]]
		write.csv(combined_data, paste0("/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output/Esse_done_results/combined_esse_runs_", taxon, ".csv"), row.names = FALSE)
	}

}


# combining the combinations with bad ess with the runs that did not converge as those are the ones I need to rerun in some way
rerun_list <- c(names(combinations_with_bad_ess), did_not_converge_null_combinations)

rerun_list
