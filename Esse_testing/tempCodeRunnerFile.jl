# Running the ESSE model
time_infer = @elapsed Tapestree.esse(tree_file, out_file, 2,
		envdata_file = envdata_file,
		states_file  = states_file, 
		out_states   = out_states,
		cov_mod      = ("s",),
		ncch         = 3,
		parallel     = true,
		niter        = 10_000,
		nthin        = 100,
		dt           = 0.8,
		nburn        = 1_000, 
		mc           = "mh",
		node_ps      = (true, 100))