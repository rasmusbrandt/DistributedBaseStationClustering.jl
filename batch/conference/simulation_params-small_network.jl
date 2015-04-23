simulation_params["I"] = 8
simulation_params["Kc"] = 1
simulation_params["N"] = 2
simulation_params["M"] = 2
simulation_params["d"] = 1

l = sqrt(simulation_params["I"]/BS_density)
simulation_params["geography_size"] = (l, l)

simulation_params["simulation_name"] = "raw-small_network"

unshift!(simulation_params["assignment_methods"], ExhaustiveSearchClustering)
