using Tapestree, DelimitedFiles, PANDA, DataFrames


# set working directory
cd("/home/au543206/Documents/TroRaiMo/Esse_testing")

# Loading tree
tree = read_newick("/home/au543206/Documents/TroRaiMo/Esse_testing/pruned_tree_order_Laurales_GBMB.tre")

# Load the array of floats from the sampling_freq file
sampling_freq_array = readdlm("Laurales_sampling_fraction.txt")
sampling_freq_array = sampling_freq_array[2:end,:] # Remove the first row with column names

sp  = convert(Vector{String},  sampling_freq_array[:,1])
sp  = replace.(sp, " " => "_")
rho = convert(Vector{Float64}, sampling_freq_array[:,2])
rho = Dict(sp[i] => rho[i] for i in 1:lastindex(sp))



name_list_tips = tip_labels(tree)

# Replacing "_" with " " in each entry in name_list_tips
name_list_tips = replace.(name_list_tips, "_" => " ")

# Add a column to name_list_tips with the numbers 1 to the number of tips in the tree.
name_list_tips = hcat(name_list_tips, 1:n_tips(tree))

# Converting to DataFrames
sampling_freq = sampling_freq[2:end,:] # Remove the first row with column names
sampling_freq = DataFrame(sampling_freq, [:species, :frequency])
name_list_tips = DataFrame(name_list_tips, [:species, :nr_in_tree])

#println(sampling_freq)
#println(name_list_tips)

# Join the sampling_freq array with the name_list_tips array by matching the first column in each array.
sampling_freq_joined = innerjoin(sampling_freq, name_list_tips, on = "species")

# pull out the frequency column as an array
sampling_freq_array = sampling_freq_joined[!, :frequency]

sampling_freq_joined


sp  = convert(Vector{String},  sampling_freq_array[:,1])
sp  = replace.(sp, " " => "_")

sp = sampling_freq_array[:,1]

rho = convert(Vector{Float64}, sampling_freq_array) # sampling frequency per species
rho = Dict(sp[i] => rho[i] for i in 1:lastindex(sp))

tr = tree

r, tv = insane_cbd(tree,
                   nburn   = 1_000,
                   niter   = 100_000,
                   nthin   = 1_000,
                   nflush  = 1_000,
                   ofile   = "/home/au543206/Documents/TroRaiMo/Esse_testing_cbd",
                   tρ      = rho,
                   mxthf   = 0.05)



iread("../Esse_testing_cbd.txt")

