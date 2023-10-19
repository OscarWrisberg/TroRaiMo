# This script adds 0's and 1's to the end of a file to represent the state of each species
# It does so at random in order for me to test how ESSE handles different amounts of states on the GBMB tree.

import random
import argparse

# Command line argument which is the file name
parser = argparse.ArgumentParser()

# Define command line argument
parser.add_argument("input_file", help="The name of the file to add a state to")
parser.add_argument("states", help="the number of states to be added to the file")
parser.add_argument("output_file", help="The name of the file to write the new states to")

# Parse the command-line arguments
args = parser.parse_args()

# open the file for reading
with open(args.input_file, 'r') as input_file:
	# read the file into a list
	lines = input_file.readlines()

# Generate a vector with 1's and 0's at random the same number of times as specified by the states argument in the command line
states = [random.randint(0,1) for i in range(int(args.states))]
print("Generated states:", states)

# Add the states to the list of species
lines.append(' '.join([str(state) for state in states]))
print("New lines:", lines)

# Open the output file for writing 
with open(args.output_file, 'w') as output_file:
	# Write the list of species to the output file
	output_file.writelines(lines)

