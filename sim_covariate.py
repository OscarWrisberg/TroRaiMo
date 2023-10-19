# This script creates a file with the covariates for the number a specific number of regions
# It does so at random in order for me to test how ESSE handles different amounts of states on the GBMB tree.

import numpy as np
import argparse

# Command line argument which is the file name
parser = argparse.ArgumentParser()

# Define command line argument
parser.add_argument("input_file", help="The name of the file to add a state to")
parser.add_argument("states", help="the number of states to be added to the file")
parser.add_argument("output_file", help="The name of the file to write the new states to")


# Parse the command-line arguments
args = parser.parse_args()

# I likely want to sim 1000 covariates for each region.
# I dont want the covariates for each region to be completely random, but I want them to be random enough that they are not correlated.
# The output file should be a txt file with the covariates for each region in a column, while the first column is the time increasing from 0 to 10.

# Define the size of the sample
sample_size = 1000

# Create a file for writing the covariates to
with open(args.output_file, 'w') as output_file:
	# Write the header line
	output_file.write('time\t')
	# Write the names of the covariates
	for i in range(int(args.states)):
		output_file.write('covariate' + str(i) + '\t')
	# Write a new line
	output_file.write('\n')
	# Write the time column
	for i in range(sample_size):
		output_file.write(str(i) + '\t')
		# Generate random numbers from a uniform distribution between 0 and 5
		random_numbers = np.random.uniform(0, 5, size=int(args.states))
		# Write the random numbers to the file
		for random_number in random_numbers:
			output_file.write(str(random_number) + '\t')
		# Write a new line
		output_file.write('\n')

