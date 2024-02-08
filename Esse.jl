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
using PANDA

# Srun file locations
# processors = 2
# tree_file = "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/orders/pruned_tree_order_Arecales_GBMB.tre"
# tip_states_file = "/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data/Arecales_states_0.2.txt"
# paleo_clim_file = "/home/au543206/GenomeDK/Trf_models/TroRaiMo/paleoclim_area.txt"
# out_states_file = "Esse_states_Arecales_0.4.txt"
# out_file = "Esse_output_"+orders[i]+".jld2"
# hidden_states = 2

# Take 6 command line arguments
processors = parse(Int, ARGS[1])
tree_file = ARGS[2]
tip_states_file = ARGS[3]
paleo_clim_file = ARGS[4]
out_states_file = ARGS[5]
out_file = ARGS[6]
save_file = ARGS[7]
hidden_states = parse(Int, ARGS[8])

# The first argument is the number of processors available to the Tapestree
# Load Tapestree and Distributed
# using Distributed
# addprocs(processors)
# @everywhere using Tapestree

# Load a tree file
tree = joinpath(tree_file)

states = joinpath(tip_states_file)

# Load the paleoenvironmental data
paleo_data = joinpath(paleo_clim_file)

# Setting the out directory for states
out_states = joinpath(out_states_file)

# Setting the outdir for the MCMC data
out_file = joinpath(out_file)
out_states = joinpath(out_states_file)


# Running the ESSE model
Tapestree.esse(tree, out_file, hidden_states,
 	envdata_file = paleo_data, 
	states_file = states,
	out_states = out_states,
	cov_mod = ("s",),
	parallel = true,
	mc = "mh",
	ncch= 3
	)


# Save the output in a JLD2 file
@save save_file out_file







