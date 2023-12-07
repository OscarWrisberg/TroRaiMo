# I want to run ClaDs on 4 different trees where 
# Arecales
# Zingiberales
# Vitales
# Geraniales

# Maybe include Fagales as there are only 15 species without dist data and it is a good temperate tree group
# Fagales (Maybe)
# Fabales is only missing 116


# path to the 4 different trees
using Pkg
Pkg.add("PANDA")
using PANDA

# arecales_tree = load_tree("/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/pruned_tree__order_Arecales_GBMB.tre")
# zingiberales_tree = load_tree("/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/pruned_tree__order_Zingiberales_GBMB.tre")
# vitales_tree = load_tree("/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/pruned_tree__order_Vitales_GBMB.tre")
# geraniales_tree = load_tree("/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/twice_pruned_tree_Geraniales_GBMB.tre")

println(ARGS[1])
println(ARGS[2])

path_to_tree = ARGS[1]
sampling_freq = ARGS[2]

tree = load_tree(path_to_tree)
output_name = path_to_tree*"_output"

output = infer_ClaDS(tree, print_state = 100, f = sampling_freq)

# output_arecales = infer_ClaDS(arecales_tree, print_state = 100, f = 0.2787973)
# output_zingiberales = infer_ClaDS(zingiberales_tree, print_state = 100, f = 0.229433)
# output_vitales = infer_ClaDS(vitales_tree, print_state = 100, f = 0.2864078)
# output_geraniales = infer_ClaDS(geraniales_tree, print_state = 100, f = 0.4295775)

using JLD2
@save "/home/owrisberg/Trf_models/workflow/03_distribution_data/"+output_name output