using Pkg


Pkg.build("DiffEqNoiseProcess")
Pkg.build("PyCall")

Pkg.add("Tapestree")


using Tapestree

tree_file = joinpath(dirname(pathof(Tapestree), "..","data","tree_50.tre"))