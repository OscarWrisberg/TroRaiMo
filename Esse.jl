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
processors = 3
tree_file = "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/orders/pruned_tree_order_Cornales_GBMB.tre"
tip_states_file = "/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data/Cornales_states_0.1.txt"
paleo_clim_file = "/home/au543206/GenomeDK/Trf_models/TroRaiMo/paleoclim_area.txt"
out_states_file = "Esse_states_Cornales_0.1.txt"
out_file = "Esse_output_Cornales_0.1.jld2"
hidden_states = 2
output_folder = "/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output/"

# path_to_tree = "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/subset_of_orders/family_phylo_Resedaceae.tre"
# sampling_freq_file = "/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data/Resedaceae_sampling_fraction.txt"


# Take command line arguments
processors = parse(Int, ARGS[1])
tree_file = ARGS[2]
tip_states_file = ARGS[3]
paleo_clim_file = ARGS[4]
out_states_file = ARGS[5]
out_file = ARGS[6]
save_file = ARGS[7]
hidden_states = parse(Int, ARGS[8])
output_folder = ARGS[9]


addprocs(processors)
@everywhere using Tapestree

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



# Running the ESSE model in parallel 
time_infer = @elapsed Tapestree.esse(tree, # Full path to the tree
	out_file, # Full path to write the MCMC output
	hidden_states, # Number of hidden states
 	envdata_file = paleo_data, # Data from koppen biomes
	states_file = states, # Data for the tip states of the species in the tree
	out_states = out_states, # The out states file where the states are saved
	cov_mod = ("s",), # s specifies that only speciation is affected by the covariate 
	parallel = true,
	mc = "mh", # Metropolis-Hastings
	ncch = processors, # number of chains
	niter = 5_000, # Number of iterations
	nthin = 100, # Frequency at which to record the state
	dt = 0.8, # Temperature for the annealing of the chains
	nburn = 1_000, # Number of iterations to use in the burn in
	)

# Print time to infer in hours
println("Time to infer: $time_infer seconds")


# Measure the time to save the output
time_save_output = @elapsed @save output_name output
println("Time to save output: $time_save_output seconds")


# Save the output in a JLD2 file
@save save_file out_file

