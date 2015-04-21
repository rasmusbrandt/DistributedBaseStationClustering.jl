simulation_params["I"] = 8
simulation_params["Kc"] = 1
simulation_params["N"] = 2
simulation_params["M"] = 2
simulation_params["d"] = 1

simulation_params["simulation_name"] = "raw-small_network"

unshift!(simulation_params["assignment_methods"], ExhaustiveSearchClustering)
unshift!(simulation_params["assignment_methods"], Chen2014_ExhaustiveSearch)
