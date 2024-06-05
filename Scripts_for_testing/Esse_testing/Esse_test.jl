# The files for the tutorial are not present in the conda download of Tapestree
# I had to download them manually from the Github

# Code for getting ressources in terminal
#srun --account Trf_models --mem 32g --pty bash

#julia


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
using JLD2

#########################################################################################
# Testing ESSE using the built in data

#Loading the tree
tree_file = joinpath("/home/au543206/esse_data/tree_50.tre")
println(tree_file)

# Loading the states
states_file = joinpath("/home/au543206/esse_data/st2_data.txt")
println(states_file)

# Loading the co-variate data
envdata_file = joinpath("/home/au543206/esse_data/env_data_2.txt")
println(envdata_file)

# Specifying out folder
out_file = joinpath("/home/au543206/esse_data/file_out_gwf.txt")
out_states = joinpath("/home/au543206/esse_data/states_out_gwf.txt")

# Running the ESSE model without parallelization
time_infer = @elapsed Tapestree.esse(tree_file,
	out_file,
	2,
	envdata_file = envdata_file, 
    states_file = states_file,
    out_states = out_states,
    cov_mod = ("s",)
	)

println("Time to run Esse: $time_infer seconds") # 3978.85 seconds or 66.3 minutes


# Running the ESSE model with parallelization
addprocs(3)
@everywhere using Tapestree

time_infer = @elapsed Tapestree.esse(tree_file,
	out_file,
	2,
	envdata_file = envdata_file, 
    states_file = states_file,
    out_states = out_states,
    cov_mod = ("s",),
	parallel = true,
	mc = "mh", # Metropolis-Hastings
	ncch = 3, # number of chains
	niter = 5_000, # Number of iterations
	nthin = 100, # Frequency at which to record the state
	dt = 0.8, # Temperature for the annealing of the chains
	nburn = 1_000)

println("Time to run Esse: $time_infer seconds")

rmprocs(3)
@everywhere using Tapestree

#########################################################################################

time_infer = @elapsed Tapestree.esse(tree_file, out_file, 2,
		envdata_file = envdata_file,
		states_file  = states_file, 
		out_states   = out_states,
		cov_mod      = ("s",),
		ncch         = 3,
		parallel     = true,
		niter        = 10_000,
		nthin        = 100,
		dt           = 0.8,
		nburn        = 1_000, 
		mc           = "mh",
		node_ps      = (true, 100))

println("Time to run Esse: $time_infer seconds") # 1419.62  seconds or 23.6 minutes


#########################################################################################
# Testing ESSE using some of my data

using Pkg
# Check if Tapestree and Distributed are installed
if !haskey(Pkg.installed(), "Tapestree") || !haskey(Pkg.installed(), "Distributed") || !haskey(Pkg.installed(), "DataFrames") || !haskey(Pkg.installed(), "DelimitedFiles") || !haskey(Pkg.installed(), "PANDA")
	# Install Tapestree and Distributed
	Pkg.add(["Tapestree", "Distributed", "DataFrames", "DelimitedFiles", "PANDA"])
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

# Take command line arguments
cd("/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/orders/")
tree_file = "pruned_tree_order_Santalales_GBMB.tre"
tip_states_file = "/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data/Santalales_states_0.2.txt"
paleo_clim_file = "/home/au543206/GenomeDK/Trf_models/data/paleoclim_area.txt"
out_states_file = "Esse_states_Santalales_0.2"
out_file = "Esse_output_Santalales_hidden_states_0.2"
save_file = "Esse_output_Santalales_0.2.jld2"
hidden_states = 2
output_folder = "/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output/"
save_file = "Esse_output_Santanales_0.2.jld2"


# The first argument is the number of processors available to the Tapestree
# Load Tapestree and Distributed
# using Distributed
# addprocs(processors)
# @everywhere using Tapestree

# Looking at the data and the tree
tree_file = joinpath(tree_file)
# Load tree file
tree = load_tree(tree_file)
name_list_tips = tip_labels(tree)

# Load the states
states_file = joinpath(tip_states_file)
states = readdlm(states_file)


# Load a tree file
tree = joinpath(tree_file)

states = joinpath(tip_states_file)

# Load the paleoenvironmental data
paleo_data = joinpath(paleo_clim_file)

# Setting the out directory for states
out_states = joinpath(output_folder,out_states_file)

# Setting the outdir for the MCMC data
out_file = joinpath(output_folder, out_file)
out_states = joinpath(output_folder, out_states_file)

# Running Esse with the following parameters
println("Running ESSE with the following parameters:")
println("Tree file: $tree")
println("out_file: $out_file")
println("Hidden states: $hidden_states")
println("Paleo data file: $paleo_data")
println("Tip states file: $states")
println("Out states file: $out_states")



# Running the ESSE model
time_infer = @elapsed Tapestree.esse(tree_file, out_file, 2,
		envdata_file = envdata_file,
		states_file  = states_file, 
		out_states   = out_states,
		cov_mod      = ("s",),
		ncch         = 3,
		parallel     = true,
		niter        = 10_000,
		nthin        = 100,
		dt           = 0.8,
		nburn        = 1_000, 
		mc           = "mh",
		node_ps      = (true, 100)) # 84721.663403089 seconds or 23.5 hours


# Save the output in a JLD2 file
@save save_file out_file

