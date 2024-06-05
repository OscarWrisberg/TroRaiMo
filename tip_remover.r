# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments are provided
if (length(args) != 3) {
	stop("Please provide the tree file and the distribution data file as command line arguments.")
}

# Extract the file paths from the command line arguments
tree_file <- args[1]
distribution_file <- args[2]
output_file <- args[3]

# Load required packages
if (!require(ape)) {
	install.packages("ape")
	library(ape)
}

if (!require(data.table)) {
	install.packages("data.table")
	library(data.table)
}

# Read the tree file
tree <- read.tree(tree_file)

# Read the distribution data file
distribution_data <- fread(distribution_file)

# Get the tips from the tree
tree_tips <- tree$tip.label

# Get the tips from the distribution data
distribution_tips <- distribution_data$V1

# Find tips in the tree that are not in the distribution data
tips_to_remove <- setdiff(tree_tips, distribution_tips)

# Printing the tips which are not in the distribution data
cat("Tips not in the distribution data:", tips_to_remove, "\n")

# Remove the tips from the tree
tree_pruned <- drop.tip(tree, tips_to_remove)

# Save the pruned tree to a new file
write.tree(tree_pruned, file = output_file)

# Print a message indicating the successful removal of tips
cat("Tips removed from the tree and the pruned tree is saved as", output_file, "\n")