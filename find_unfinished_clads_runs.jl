#########################################################################################################################
################################################## Loading packages  ####################################################
#########################################################################################################################
# Load the required packages
using Pkg
# Check if Tapestree and Distributed are installed
if !haskey(Pkg.installed(), "DataFrames") || !haskey(Pkg.installed(), "FileIO") || !haskey(Pkg.installed(), "JLD2") || !haskey(Pkg.installed(), "GlobalScope") || !haskey(Pkg.installed(), "Printf")
	# Install Tapestree and Distributed
	Pkg.add(["DataFrames", "FileIO", "JLD2", "Glob", "Printf"])
end

# Load Tapestree and Distributed
using DataFrames
using FileIO
using JLD2
using Glob
using Printf
using PANDA

# Setting output folder
output_folder = "/home/au543206/GenomeDK/Trf_models/workflow/05_figures"

########################################################################################################################
############################################## Checking for duplicates  ################################################
########################################################################################################################

# Clads on orders
folderpath_1 = "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/orders"
file_list_1 = glob("Clads_output_*.jld2", folderpath_1)

# Clads on suborders
folderpath_2 = "/home/au543206/GenomeDK/Trf_models/workflow/02_adding_orders/pruning/subset_of_orders"
file_list_2 = glob("Clads_output_*.jld2", folderpath_2)

# Join the file lists
file_list = vcat(file_list_1, file_list_2)

println(file_list)  # length should be 264 as of 30/08/2024

########################################################################################################################
############################################## Creating the dataframe  ################################################
########################################################################################################################
Unfinished_runs = String[]

# Loop through each file in the folder
for (i, file) in enumerate(file_list)

    # Keeping track of the progress
    println("Processing file $i of $(length(file_list))")

    # Extracting the order name from the file name
    order_name = replace(basename(file), r"Clads_output_(.*)\.jld2" => s"\1")

    # Load the JLD2 file
    CladsOutput = load(file)

    if haskey(CladsOutput, "reason_for_stop")
        println("The reason for stop is in the file")
        if CladsOutput["output"].reason_for_stop == "end_time"
            println("The reason for stop is because the run ran out of time")
            push!(Unfinished_runs, order_name)
        end
    else
        continue
    end

    # Clean up the CladsOutput object
    empty!(CladsOutput)
end

# Printing the list of unfinished runs
println(Unfinished_runs)


##########################################
# Find the Cucurbitaceae run in file list
Cucurbitaceae_run = findfirst(x -> occursin("Cucurbitaceae", x), file_list)


 CladsOutput = load(file_list[Cucurbitaceae_run])
 println(keys(CladsOutput))
 println(CladsOutput["tree"]["tip.label"])

# Print the Cladsoutput output
println(CladsOutput["output"])

#
CladsOutput["output"].reason_for_stop