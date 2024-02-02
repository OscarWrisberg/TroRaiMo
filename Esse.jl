using Pkg
# Check if Tapestree and Distributed are installed
if !haskey(Pkg.installed(), "Tapestree") || !haskey(Pkg.installed(), "Distributed") || !haskey(Pkg.installed(), "DataFrames") || !haskey(Pkg.installed(), "DelimitedFiles")
	# Install Tapestree and Distributed
	Pkg.add(["Tapestree", "Distributed", "DataFrames", "DelimitedFiles"])
end

# Load Tapestree and Distributed
using Tapestree
using Distributed
using DataFrames
using DelimitedFiles


# Srun file locations
# processors = 2
# tree_file = "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/orders/pruned_tree_order_Arecales_GBMB.tre"
# tip_states_file = "/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data/Arecales_distribution_data.txt"
# paleo_clim_file = "/home/au543206/GenomeDK/Trf_models/TroRaiMo/paleoclim_area.txt"
# out_states_file = "Esse_states_Laurales_0.4.txt"
# out_file = "whatever2"
# hidden_states = 2
# percentage_for_present = 0.2

# Take 6 command line arguments
processors = parse(Int, ARGS[1])
tree_file = ARGS[2]
tip_states_file = ARGS[3]
paleo_clim_file = ARGS[4]
out_states_file = ARGS[5]
out_file = ARGS[6]
hidden_states = parse(Int, ARGS[7])

# The first argument is the number of processors available to the Tapestree
# Load Tapestree and Distributed
using Distributed
addprocs(processors)
@everywhere using Tapestree

# Load a tree file
tree = joinpath(tree_file)


states = joinpath()

# Load the paleoenvironmental data
paleo_data = joinpath(paleo_clim_file)
paleo_data


# Setting the out directory for states
out_states = joinpath(out_states_file)

# Setting the outdir for the MCMC data
out = joinpath(out_file)

# If there are any tips in the tree which are not in the distribution data, drop them.

esse(tree, out, hidden_states, envdata_file = paleo_data, 
  states_file = states, out_states = out_states, cov_mod = ("s",), parallel = true)


# Saving the output in a JLD2 file.
@save output_name output

# Savint the output in R
r_output = replace(output_name, r"\.jld2" => ".R")





