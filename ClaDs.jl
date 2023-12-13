# Loading packages
using PANDA
using JLD2

# Measure the time to load Pkg
#time_load_pkg = @elapsed using Pkg

# Measure the time to load PANDA
time_load_panda = @elapsed using PANDA

# Measure the time to load JLD2
time_load_jld2 = @elapsed using JLD2

# Print the time to load the packages
println("Time to load PANDA: $time_load_panda seconds")
println("Time to load JLD2: $time_load_jld2 seconds")

# Prepare some paths test script locally
# path_to_tree = "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/pruned_tree__order_Zingiberales_GBMB.tre"
# output_name = "test_output.jld2"

# Fetching arguments
path_to_tree = ARGS[1]
sampling_freq = parse(Float64, ARGS[2])
output_name = ARGS[3]


# Measure the time to load the tree
time_load_tree = @elapsed tree = load_tree(path_to_tree)
println("Time to load the tree: $time_load_tree seconds")

# Measure the time to run infer_ClaDS
time_infer = @elapsed output = infer_ClaDS(tree, print_state = 100, f = sampling_freq)
println("Time to run infer_ClaDS: $time_infer seconds")

# Measure the time to save the output
time_save_output = @elapsed @save output_name output
println("Time to save output: $time_save_output seconds")
