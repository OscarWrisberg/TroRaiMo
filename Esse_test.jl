# The files for the tutorial are not present in the conda download of Tapestree
# I had to download them manually from the Github

# Code for getting ressources in terminal
#srun --account Trf_models --mem 32g --pty bash

#julia

# Loading Tapestree
using Tapestree

#Loading the tree
tree_file = joinpath(dirname(pathof(Tapestree)), "..", "data", "tree_50.tre")
println(tree_file)

# Loading the states
states_file = joinpath(dirname(pathof(Tapestree)), "..", "data", "st2_data.txt")
println(states_file)

# Loading the co-variate data
envdata_file = joinpath(dirname(pathof(Tapestree)), "..", "data", "env_data_2.txt")
println(envdata_file)

# Specifying out folder
out_file = *(homedir(),"/Trf_models/Esse_test/file_out_gwf.txt")
out_states = *(homedir(),"/Trf_models/Esse_test/states_out_gwf.txt")

esse(tree_file, out_file, 2, envdata_file = envdata_file, 
  states_file = states_file, out_states = out_states, cov_mod = ("s",))
