using Pkg
# Check if Tapestree and Distributed are installed
if !haskey(Pkg.installed(), "PANDA") || !haskey(Pkg.installed(), "JLD2") || !haskey(Pkg.installed(), "DataFrames") || !haskey(Pkg.installed(), "DelimitedFiles")
	# Install Tapestree and Distributed
	Pkg.add(["PANDA", "JLD2", "DataFrames", "DelimitedFiles"])
end

# Measure the time to load Pkg
#time_load_pkg = @elapsed using Pkg

# Measure the time to load PANDA
time_load_panda = @elapsed using PANDA

# Measure the time to load JLD2
time_load_jld2 = @elapsed using JLD2

# Measure the time to load DataFrames
time_load_dataframes = @elapsed using DataFrames

# Measure the time to load DelimitedFiles
time_load_delimitedfiles = @elapsed using DelimitedFiles

# Print the time to load the packages
println("Time to load PANDA: $time_load_panda seconds")
println("Time to load JLD2: $time_load_jld2 seconds")
println("Time to load DataFrames: $time_load_dataframes seconds")
println("Time to load DelimitedFiles: $time_load_delimitedfiles seconds")


# Prepare some paths test script locally
path_to_tree = "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/orders/pruned_tree_order_Zingiberales_GBMB.tre"
sampling_freq_file = "/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data/Zingiberales_sampling_fraction.txt"

path_to_tree = "/home/owrisberg/Trf_models/workflow/02_adding_orders/pruning/subset_of_orders/family_phylo_Lauraceae.tre"
sampling_freq_file = "/home/owrisberg/Trf_models/workflow/03_distribution_data/Laurales_sampling_fraction.txt"


output_name = "test_output.jld2"


# Fetching arguments
path_to_tree = ARGS[1]
sampling_freq_file = ARGS[2]
output_name = ARGS[3]

# Measure the time to load the tree
time_load_tree = @elapsed tree = load_tree(path_to_tree)
println("Time to load the tree: $time_load_tree seconds")

# Load the array of floats from the sampling_freq file
sampling_freq = readdlm(sampling_freq_file)
name_list_tips = tip_labels(tree)

# Replacing "_" with " " in each entry in name_list_tips
name_list_tips = replace.(name_list_tips, "_" => " ")

# Add a column to name_list_tips with the numbers 1 to the number of tips in the tree.
name_list_tips = hcat(name_list_tips, 1:n_tips(tree))

# Converting to DataFrames
sampling_freq = sampling_freq[2:end,:] # Remove the first row with column names
sampling_freq = DataFrame(sampling_freq, [:species, :frequency])
name_list_tips = DataFrame(name_list_tips, [:species, :nr_in_tree])

println(sampling_freq)
#println(name_list_tips)

# Join the sampling_freq array with the name_list_tips array by matching the first column in each array.
sampling_freq_joined = innerjoin(sampling_freq, name_list_tips, on = "species")

# pull out the frequency column as an array
sampling_freq_array = sampling_freq_joined[!, :frequency]

# Measure the time to run infer_ClaDS
time_infer = @elapsed output = infer_ClaDS(tree,
	print_state = 100,
	f = sampling_freq_array,)
	println("Time to run infer_ClaDS: $time_infer seconds")

# Measure the time to save the output
time_save_output = @elapsed @save output_name output
println("Time to save output: $time_save_output seconds")

r_output = replace(output_name, r"\.jld2" => ".Rdata")

# Also save the output to an R file
save_ClaDS_in_R(output,  r_output)
