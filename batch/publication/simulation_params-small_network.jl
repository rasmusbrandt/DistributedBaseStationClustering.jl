simulation_params["I"] = 8
simulation_params["Kc"] = 1
simulation_params["N"] = 2
simulation_params["M"] = 2
simulation_params["d"] = 1

unshift!(simulation_params["assignment_methods"], BranchAndBoundClustering)
unshift!(simulation_params["assignment_methods"], Chen2014_ExhaustiveSearch)
