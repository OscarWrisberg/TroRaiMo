using Pkg

using Tapestree

tree_file = joinpath("/home/au543206/Tapestree.jl/data/tree_50.tre")

states_file = joinpath("/home/au543206/Tapestree.jl/data/st2_data.txt")

env_file = joinpath("/home/au543206/Tapestree.jl/data/env_data_2.txt")

out_file = *(homedir(),"/Test_1_chain")

out_states = *(homedir(),"/Test_1_chain.logs")

cpar = ("q_01 = q_10",
        "loss_A_0 = mu_A_0",
        "loss_B_0 = mu_B_0",
        "loss_A_1 = mu_A_1",
        "loss_B_1 = mu_B_1")

esse(tree_file,
		out_file,
		2,
		envdata_file = env_file,
		states_file  = states_file, 
		out_states   = out_states,
		constraints  = cpar,
		cov_mod      = ("s",),
		niter        = 5_000,
		nthin        = 100,
		dt           = 0.8,
		nburn        = 1_000, 
		mc           = "mh",
		node_ps      = (true, 100))


using Distributed
addprocs(3)
@everywhere using Tapestree

out_file = *(homedir(),"/Test_3_chain")

out_states = *(homedir(),"/Test_3_chain.logs")


esse(tree_file, out_file, 2,
     envdata_file = env_file,
     states_file  = states_file, 
     out_states   = out_states,
     constraints  = cpar,
     cov_mod      = ("s",),
     ncch         = 3,
     parallel     = true,
     niter        = 5_000,
     nthin        = 100,
     dt           = 0.8,
     nburn        = 1_000, 
     mc           = "mh",
     node_ps      = (true, 100))