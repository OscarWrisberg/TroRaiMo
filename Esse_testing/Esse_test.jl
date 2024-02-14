# The files for the tutorial are not present in the conda download of Tapestree
# I had to download them manually from the Github

# Code for getting ressources in terminal
#srun --account Trf_models --mem 32g --pty bash

#julia

#########################################################################################
# Testing ESSE using the built in data
# Loading Tapestree
using Tapestree

#Loading the tree
tree_file = joinpath(dirname(pathof(Tapestree)), "..", "data", "tree_50.tre")
println(tree_file)

# Loading the states
states_file = joinpath(dirname(pathof(Tapestree)), "..", "data", "st2_data.txt")
println(states_file)

# Loading the co-variate data
envdata_file = joinpath(dirname(pathof(Tapestree)), "..", "data", "env_data_2.txt")
println(envdata_file)

# Specifying out folder
out_file = *(homedir(),"/Trf_models/Esse_test/file_out_gwf.txt")
out_states = *(homedir(),"/Trf_models/Esse_test/states_out_gwf.txt")

esse(tree_file, out_file, 2, envdata_file = envdata_file, 
  states_file = states_file, out_states = out_states, cov_mod = ("s",))



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
cd("/home/owrisberg/Trf_models/workflow/02_adding_orders/pruning/orders/")
processors = 1
tree_file = "pruned_tree_order_Santalales_GBMB.tre"
tip_states_file = "/home/owrisberg/Trf_models/workflow/03_distribution_data/Santalales_states_0.2.txt"
paleo_clim_file = "/home/owrisberg/Trf_models/data/paleoclim_area.txt"
out_states_file = "Esse_states_Santalales_0.2"
out_file = "Esse_output_Santalales_hidden_states_0.2"
save_file = "Esse_output_Santalales_0.2.jld2"
hidden_states = 2
output_folder = "/home/owrisberg/Trf_models/workflow/04_results/Esse_output/"

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
Tapestree.esse(tree, out_file, hidden_states
 	envdata_file = paleo_data, 
	states_file = states,
	out_states = out_states,
	cov_mod = ("s",)
	)


# Save the output in a JLD2 file
@save save_file out_file


