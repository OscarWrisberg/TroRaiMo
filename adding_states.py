# # This script adds 0's and 1's to the end of a file to represent the state of each species
# # It does so at random in order for me to test how ESSE handles different amounts of states on the GBMB tree.

# import random
# import argparse

# # Command line argument which is the file name
# parser = argparse.ArgumentParser()

# # Define command line argument
# parser.add_argument("input_file", help="The name of the file to add a state to")
# parser.add_argument("states", help="the number of states to be added to the file")
# parser.add_argument("output_file", help="The name of the file to write the new states to")

# # Parse the command-line arguments
# args = parser.parse_args()

# # open the file for reading
# with open(args.input_file, 'r') as input_file:
# 	# read the file into a list
# 	lines = input_file.readlines()

# # Generate a vector with 1's and 0's at random the same number of times as specified by the states argument in the command line and the number of species
# for j in lines:
# 	states = [random.randint(0,1) for i in range(int(args.states))]
# 	print("Generated states:", states)
# 	# Add the states to the list of species
# 	j.append(' '.join([str(state) for state in states]))
# 	print("Species with states:", j)



# # Open the output file for writing 
# with open(args.output_file, 'w') as output_file:
# 	# Write the list of species to the output file
# 	output_file.writelines(lines)

# This script adds 0's and 1's to each species across all columns to represent the state of each species.

import random
import argparse

# Command line argument which is the file name
parser = argparse.ArgumentParser()

# Define command line argument
parser.add_argument("input_file", help="The name of the file to add states to")
parser.add_argument("states", type=int, help="the number of states to add to each species")
parser.add_argument("output_file", help="The name of the file to write the new states to")

# Parse the command-line arguments
args = parser.parse_args()

# open the file for reading
with open(args.input_file, 'r') as input_file:
    # read the file into a list
    lines = input_file.readlines()

# Determine the number of species
num_species = len(lines)

# Generate random 0's and 1's for each species and add them to each column
output_lines = []
for _ in range(num_species):
    # Generate random 0's and 1's for the specified number of states
    states = [str(random.randint(0, 1)) for _ in range(args.states)]
    
    # Join the states and append to the output
    output_lines.append(' '.join(states))

# Open the output file for writing 
with open(args.output_file, 'w') as output_file:
    # Write the list of species with added states to the output file
    output_file.writelines('\n'.join(output_lines))


