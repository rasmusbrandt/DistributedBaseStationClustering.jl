simulation_params["I"] = 16
simulation_params["Kc"] = 1
simulation_params["N"] = 2
simulation_params["M"] = 4
simulation_params["d"] = 1

l = sqrt(simulation_params["I"]/BS_density)
simulation_params["geography_size"] = (l, l)
