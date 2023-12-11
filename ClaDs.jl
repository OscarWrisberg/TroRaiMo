# Loading packages
using Pkg
using PANDA
using JLD2

# Measure the time to load Pkg
time_load_pkg = @elapsed using Pkg

# Measure the time to load PANDA
time_load_panda, _ = @time using PANDA

# Measure the time to load JLD2
time_load_jld2, _ = @time using JLD2

# Print the time to load the packages
println("Time to load Pkg: $time_load_pkg seconds")
println("Time to load PANDA: $time_load_panda seconds")
println("Time to load JLD2: $time_load_jld2 seconds")

# Fetching arguments
path_to_tree = ARGS[1]
sampling_freq = parse(Float64, ARGS[2])

# Measure the time to load the tree
time_load_tree, _ = @time tree = load_tree(path_to_tree)
println("Time to load the tree: $time_load_tree seconds")

# Measure the time to run infer_ClaDS
time_infer, _ = @time output = infer_ClaDS(tree, print_state = 100, f = sampling_freq)
println("Time to run infer_ClaDS: $time_infer seconds")

# Measure the time to save the output
time_save_output, _ = @time @save ARGS[3] output
println("Time to save output: $time_save_output seconds")



###################################################################################################################################################3
# Old stuff

# arecales_tree = load_tree("/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/pruned_tree__order_Arecales_GBMB.tre")
# zingiberales_tree = load_tree("/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/pruned_tree__order_Zingiberales_GBMB.tre")
# vitales_tree = load_tree("/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/pruned_tree__order_Vitales_GBMB.tre")
# geraniales_tree = load_tree("/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/twice_pruned_tree_Geraniales_GBMB.tre")


# output_arecales = infer_ClaDS(arecales_tree, print_state = 100, f = 0.2787973)
# output_zingiberales = infer_ClaDS(zingiberales_tree, print_state = 100, f = 0.229433)
# output_vitales = infer_ClaDS(vitales_tree, print_state = 100, f = 0.2864078)
# output_geraniales = infer_ClaDS(geraniales_tree, print_state = 100, f = 0.4295775)