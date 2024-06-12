using Pkg
# Check if Pacjages are installed
if !haskey(Pkg.installed(), "JLD2") || !haskey(Pkg.installed(), "DataFrames") || !haskey(Pkg.installed(), "DelimitedFiles") || !haskey(Pkg.installed(), "Revise")
    # Install Tapestree and Distributed
    Pkg.add(["JLD2", "DataFrames", "DelimitedFiles","Revise"])
end


# Testing the the clads_output.jl file
using Revise

# Measure the time to load JLD2 
using JLD2

# Measure the time to load DataFrames
using DataFrames

# Measure the time to load DelimitedFiles
using DelimitedFiles

# Load the local version of PANDA
Pkg.activate("/home/au543206/Documents/github_reps/PANDA.jl") # ] activate "/home/au543206/Documents/github_reps/PANDA.jl"
Pkg.instantiate() # ] instantiate
using PANDA



# Prepare some paths test script locally
# path_to_tree = "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/orders/pruned_tree_order_Zingiberales_GBMB.tre"
# sampling_freq_file = "/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data/Zingiberales_sampling_fraction.txt"

# path_to_tree = "/home/owrisberg/Trf_models/workflow/02_adding_orders/pruning/subset_of_orders/family_phylo_Lauraceae.tre"
# sampling_freq_file = "/home/owrisberg/Trf_models/workflow/03_distribution_data/Laurales_sampling_fraction.txt"

path_to_tree = "/home/au543206/esse_data/tree_50.tre"


output_name = "test_output.jld2"

# Measure the time to load the tree
tree = load_tree(path_to_tree)

# run infer_ClaDS
output = infer_ClaDS(tree,
    print_state=100,
    end_tme=1)

# Print the reason for stopping
output.reason_for_stop

# Rerun ClaDs from last output
output_2 = infer_ClaDS(tree,
					former_run=output,
					print_state=100,
					end_tme=1)

# Rerun ClaDs from last output
output_3 = infer_ClaDS(tree,
					former_run=output_2,
					print_state=100,
					end_tme=10)

output_3.reason_for_stop # Goal gelman!
