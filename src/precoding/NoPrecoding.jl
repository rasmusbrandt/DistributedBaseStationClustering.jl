##########################################################################
# NoPrecoding
#
# Pseudo method that does not perform anything but returning empty results.

function NoPrecoding(channel, network)
    K = get_num_MSs(network)
    ds = get_num_streams(network); max_d = maximum(ds)
    aux_params = get_aux_precoding_params(network)

    utilities = Array(Float64, K, max_d, 1)

    z1 = zeros(Float64, 1)
    z2 = zeros(Float64, K, max_d, 1)

    results = PrecodingResults()
    if aux_params["output_protocol"] == :all_iterations
        results["objective"] = z1
        results["utilities"] = utilities
        results["weighted_logdet_rates_full"] = z2
        results["weighted_logdet_rates_partial"] = z2
        results["allocated_power"] = z2
    elseif aux_params["output_protocol"] == :final_iteration
        results["objective"] = z1[end]
        results["utilities"] = utilities[:,:,end]
        results["weighted_logdet_rates_full"] = z2[:,:,end]
        results["weighted_logdet_rates_partial"] = z2[:,:,end]
        results["allocated_power"] = z2[:,:,end]
    end
    return results
end
