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


# Local testing
tip_states_file = "/home/au543206/GenomeDK/Trf_models/workflow/03_distribution_data/Arecales_distribution_data.txt"
out_file = "whatever2"
percentage_for_present = 0.2

# Take the command line arguments
println("Command line arguments")
tip_states_file = string(ARGS[1])
out_file = string(ARGS[2])
percentage_for_present = parse(Float64, ARGS[3])

# Load the tip states file
println("Loading the tip states file")
tip_states = readdlm(tip_states_file, '\t')
tip_states_header = tip_states[1,:]
tip_states = tip_states[2:end,:]

# Convert to a DataFrame with the column names
println("Converting to a DataFrame")
tip_states = DataFrame(tip_states, Symbol.(tip_states_header))
println(tip_states)

# Add columns based on conditions
println("Adding columns based on conditions")
transform!(tip_states, :proportion_in_tropical_rainforest => ByRow(x -> x > percentage_for_present ? 1 : 0) => :present_in_trf,
           :proportion_outside_tropical_rainforest => ByRow(x -> x > percentage_for_present ? 1 : 0) => :present_outside_trf)

# Select subset of columns
println("Selecting subset of columns")
tip_states_subset = tip_states[:,[:wcvp_taxon_name, :present_in_trf, :present_outside_trf]]

# Convert to matrix and remove the header
tip_states_subset = Matrix(tip_states_subset)

# Change " " to "_" in each species name_list_tips
println("Changing ' ' to '_' in each species name")
tip_states_subset[:,1] = replace.(tip_states_subset[:,1], " " => "_")

# Save as a tab seperated file
println("Saving as a tab seperated file")
writedlm(joinpath(out_file), tip_states_subset, '\t')

# Exit
println("Done")