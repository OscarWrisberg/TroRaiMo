
# Load the coda package
library(coda)
library(data.table)

# Setting the directory
esse_dir <- "/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output"

# Finding all the files which end in _0.33.log
files <- list.files(esse_dir, pattern = "_0.33.log", full.names = TRUE)

# Now we sort the files into 2 lists, one which have ESS above 200 and one which does not.
# We do this with a for loop
good_files <- c()
bad_files <- c()

for (file in files) {
  # Read the file
  bayesian_output <- read.table(file, header = TRUE)
  
  # Convert the data to a coda object
  coda_object <- as.mcmc(bayesian_output)
  
  # Calculate the ESS for each parameter
  ess <- effectiveSize(coda_object)

  # print ess
  print(ess)
  
  # Check if all ESS values are above 200
  if (ess[3] > 200 & ess[4] > 200 & ess[5] > 200) {
	good_files <- c(good_files, file)
  } else {
	bad_files <- c(bad_files, file)
  }
}

# Print the good files
print(good_files)

# Print the bad files
print(bad_files)

ess[3]
ess[4]
ess[5]


# Convert the data to a coda object
coda_object <- as.mcmc(bayesian_output)

# Calculate the ESS for each parameter
ess <- effectiveSize(coda_object)

# Print the ESS values
print(ess)