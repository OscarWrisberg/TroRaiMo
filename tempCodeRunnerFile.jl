# Measure the time to load the tree
time_load_tree = @elapsed tree = load_tree(path_to_tree)
println("Time to load the tree: $time_load_tree seconds")
