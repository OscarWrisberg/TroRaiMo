using Pkg
# Check if Tapestree and Distributed are installed
if !haskey(Pkg.installed(), "Tapestree") || !haskey(Pkg.installed(), "Distributed") || !haskey(Pkg.installed(), "DataFrames") || !haskey(Pkg.installed(), "DelimitedFiles")
	# Install Tapestree and Distributed
	Pkg.add(["Tapestree", "Distributed", "DataFrames", "DelimitedFiles","JLD2"])
end

# Load Tapestree and Distributed
using Tapestree
using Distributed
using DataFrames
using DelimitedFiles
using PANDA
using JLD2

# Srun file locations
processors = 3
tree_file = "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/subset_of_orders/Apiaceae_2_Esse_tree.tre"
tip_states_file = "/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data/Apiaceae_2_states_0.1.txt"
paleo_clim_file = "/home/au543206/GenomeDK/Trf_models/TroRaiMo/paleoclim_area.txt"
out_states_file = "Esse_states_Apiaceae_2_test_0.1.txt"
out_file = "Esse_output_Apiaceae_2_test_0.1.jld2"
hidden_states = 1
output_folder = "/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output/"
biome_sampling = "/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data/Apiaceae_2_biome_sampling_fraction.txt"

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
biome_sampling = ARGS[10]
niter_arg = parse(Int,ARGS[11])


# Print the command line arguments
println("Using ", processors, " processors\n")
println("Tree file: ", tree_file, "\n")
println("Tip states file: ", tip_states_file, "\n")
println("Paleo climate file: ", paleo_clim_file, "\n")
println("Output states file: ", out_states_file, "\n")
println("Output file: ", out_file, "\n")
println("Save file: ", save_file, "\n")
println("Hidden states: ", hidden_states, "\n")
println("Output folder: ", output_folder, "\n")
println("Biome sampling: ", biome_sampling, "\n")
println("Number of iterations is: ", niter_arg, "\n\n\n")

# Open file with biome sampling sample_fractions
sampling_freq = readdlm(biome_sampling, Float64, header = true)

# Removing the annoying x in the header
sampling_freq = sampling_freq[1]

#Converting to vector
sampling_freq = vec(sampling_freq)

# print the sampling fractions 
println("Sampling fractions per biome (Trf, Non-trf and Widespread): $sampling_freq")

# Add processors
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
println("\n\n Running ESSE with the following parameters: \n")
println("Tree file: $tree \n")
println("out_file: $out_file \n")
println("Hidden states: $hidden_states \n")
println("Paleo data file: $paleo_data \n")
println("Tip states file: $states \n")
println("The numbe of iterations is $niter_arg")
println("Out states file: $out_states \n\n")

# Running the ESSE model in parallel 
time_infer = @elapsed Tapestree.esse(tree, # Full path to the tree
	constraints  = ("loss_A_0 = mu_A_0","loss_B_0 = mu_B_0"), # Constraining local and global extinction rates to be the same 
	out_file, # Full path to write the MCMC output
	hidden_states, # Number of hidden states
 	envdata_file = paleo_data, # Data from koppen biomes through time
	states_file = states, # Data for the tip states of the species in the tree 
	out_states = out_states, # The out states file where the states are saved
	cov_mod = ("s",), # s specifies that only speciation is affected by the covariate 
	parallel = true, # run MCMC-mh3 in parallel
	mc = "mh", # Metropolis-Hastings
	ntakew = 200, # Number of iterations from Nburn to tune the window
	ncch = processors, # number of cores for the parallel run
	niter = niter_arg, # Number of iterations
	nthin = 100, # Frequency at which to record the state
	dt = 0.8, # Temperature for the annealing of the chains
	nburn = 20_000, # Number of iterations to use in the burn in
	ρ = sampling_freq # Sampling fractions for the biomes
	)

# Print time to infer in hours
println("Time to infer: $time_infer seconds")


# Measure the time to save the output
time_save_output = @elapsed @save save_file out_file
println("Time to save output: $time_save_output seconds")


# # Save the output in a JLD2 file
# @save save_file out_file

