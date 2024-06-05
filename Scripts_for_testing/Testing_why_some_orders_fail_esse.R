# This is a script to figure out why some Esse runs fail.

# This is the following error message that I get for the Apiales order:
# ┌ Warning: Pkg.installed() is deprecated
# └ @ Pkg ~/miniconda3/envs/Julia_env/share/julia/stdlib/v1.9/Pkg/src/Pkg.jl:745
# ┌ Warning: Pkg.installed() is deprecated
# └ @ Pkg ~/miniconda3/envs/Julia_env/share/julia/stdlib/v1.9/Pkg/src/Pkg.jl:745
#    Resolving package versions...
#   No Changes to `~/miniconda3/envs/Julia_env/share/julia/environments/Julia_env/Project.toml`
#   No Changes to `~/miniconda3/envs/Julia_env/share/julia/environments/Julia_env/Manifest.toml`
# ERROR: LoadError: Data file cannot be made of the right dimensions.
#  Make sure the data file has the same number of rows as tips in the tree
# Stacktrace:
#  [1] error(s::String)
#    @ Base ./error.jl:35
#  [2] read_data_esse(states_file::String, tree_file::String, envdata_file::String)
#    @ Tapestree.ESSE ~/.julia/packages/Tapestree/6XNCR/src/esse/sse_wrapper.jl:492
#  [3] esse(tree_file::String, out_file::String, h::Int64; states_file::String, envdata_file::String, cov_mod::Tuple{String}, node_ps::Tuple{Bool, Int64}, out_states::String, constraints::Tuple{String}, mvpars::Tuple{String}, niter::Int64, nthin::Int64, nburn::Int64, tune_int::Int64, nswap::Int64, ncch::Int64, parallel::Bool, dt::Float64, ntakew::Int64, winit::Float64, scale_y::Tuple{Bool, Bool}, algorithm::String, mc::String, λpriors::Float64, μpriors::Float64, gpriors::Float64, lpriors::Float64, qpriors::Float64, βpriors::Tuple{Float64, Float64}, hpriors::Float64, optimal_w::Float64, tni::Float64, obj_ar::Float64, screen_print::Int64, Eδt::Float64, ti::Float64, ρ::Vector{Float64})
#    @ Tapestree.ESSE ~/.julia/packages/Tapestree/6XNCR/src/esse/sse_wrapper.jl:105
#  [4] top-level scope
#    @ ./timing.jl:393
# in expression starting at /faststorage/project/Trf_models/TroRaiMo/Esse.jl:72
# srun: error: cn-1036: task 0: Exited with exit code 1

# Now I will try to figure out why this is happening.

# Start by loading the tree which is causing the error.
tree_file <- "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/orders/pruned_tree_order_Apiales_GBMB.tre"

# Load the tree
tree <- ape::read.tree(tree_file)

# Load the data file
data_file <- "/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data/Apiales_states_0.1.txt" 

# Load the data
data <- data.table::fread(data_file)

dim(data) #2404


length(tree$tip.label) #2407

# it seems like there a 3 tips which are in the tree but not in the data file.
# What are these tips? Lets use the function from Phytools to find out.
tips_not_in_data <- setdiff(tree$tip.label, data$V1)

tips_not_in_data

