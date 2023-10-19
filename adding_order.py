# In this script I want to look at all the tips in the GBMB tree, and find their order in the WCVP names file.
# I then want to cut the GBMB tree into subtrees for each of the orders present in the GBMB tree
# I then want to save each of these trees with the order name and GBMB as a part of the file name.
# Finally, I want to calculate the number of tips in each of these subtrees and save that to a file where I in the header write (Pruning method and number of trees) and then in each row following I will have the name of the tree file followed by the number of tips.

# Loading the required packages
import pandas as pd
import Bio.Phylo.NewickIO as NewickIO


# Command line argument which is the file name
parser = argparse.ArgumentParser()

# Define command line argument
parser.add_argument("input_file_tree", help="The name of the file to add a state to")
parser.add_argument("output_file", help="The name of the file to write the new states to")
parser.add_argument("input_file_wcvp", help = "The name of the file containing the WCVP names")

# Parse the command-line arguments
args = parser.parse_args()

# Open the file containing the WCVP names
with open(args.input_file_wcvp, 'r') as input_file_wcvp:
	# Read the file into a pandas dataframe, the delimiter is a |, and the header is the first row
	lines_wcvp = pd.read_csv(input_file_wcvp, sep = '|', header = 0)

# Open the file containing the GBMB tree using phylopandas
tree = NewickIO.read(args.input_file_tree)

# Find the tips in the tree
tips = tree.tips()

# Find the tips which are in the WCVP names file  
matching_tips = [tip for tip in tips if tip in lines_wcvp]

# Find the tips which are not in the WCVP names file
not_matching_tips = [tip for tip in tips if tip not in lines_wcvp]
print("These tips are not in the WCVP \n" , not_matching_tips)


# Loop through the matching tips and find the order of each tip by looking it up in the WCVP names file.
# Then append each species name and order pair to a pandas dataframe.
# Is this dataframe going to contain only one species as I erase the df_orders everytime I iterate to a new name in matching tips?
#
def find_order(name_list, wcvp):
	# Find the order of the name in WCVP names file
	names = []
	orders = []
	for name in name_list:
		# Find the order of the name in WCVP names file
		order = wcvp.loc[wcvp['taxon_name'] == name, 'order']
		print(name, order)
		# Save the name and order to a pandas dataframe
		names.append(name)
		orders.append(order)

	# Create a pandas dataframe with the names and orders
	df_orders = pd.DataFrame({'name': names, 'order': orders})
	return(df_orders)

tips_orders = find_order(matching_tips, lines_wcvp)

print(tips_orders)

# write tips orders to a txt file
tips_orders.to_csv('tips_orders.txt', sep = '\t', index = False)


# Find the unique orders in the tips_orders dataframe and which tips belong to each order
# Then prune the tree for each order and save the pruned tree to a file with the order name and GBMB as a part of the file name.

# Find the unique orders in the tips_orders dataframe
unique_orders = tips_orders.order.unique()
print("Unique orders", unique_orders)

# Find which tips belong to each order




