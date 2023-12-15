# Check if Tapestree and Distributed are installed
if !haskey(Pkg.installed(), "Tapestree") || !haskey(Pkg.installed(), "Distributed")
	# Install Tapestree and Distributed
	using Pkg
	Pkg.add(["Tapestree", "Distributed"])
end


# Take 6 command line arguments
processors = parse(Int, ARGS[1])
tree_file = parse(Int, ARGS[2])
tip_states_file = parse(Int, ARGS[3])
paleo_clim_file = parse(Int, ARGS[4])
out_states_file = parse(Int, ARGS[5])
out_file = parse(Int, ARGS[6])
hidden_states = parse(Int, ARGS[7])

# The first argument is the number of processors available to the Tapestree
# Load Tapestree and Distributed
using Distributed
addprocs(processors)
@everywhere using Tapestree

# Load a tree file
tree = joinpath(tree_file)

# Load the tip states
tip_states = joinpath(tip_states_file)

# Load the paleoenvironmental data
paleo_data = joinpath(paleo_clim_file)

# Setting the out directory for states
out_states = joinpath(out_states_file)

# Setting the outdir for the MCMC data
out = joinpath(out_file)


esse(tree, out, hidden_states, envdata_file = paleo_data, 
  states_file = tip_states, out_states = out_states, cov_mod = ("s",))

