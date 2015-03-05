function partition_to_cluster_assignment_matrix(network, partition)
    I = get_no_BSs(network); K = get_no_MSs(network)
    assignment = get_assignment(network)

    cluster_assignment_matrix = zeros(Int, K, I)

    for block in partition.blocks
        for i in block.elements
            for j in block.elements; for l in served_MS_ids(j, assignment)
                cluster_assignment_matrix[l,i] = 1
            end; end
        end
    end

    return cluster_assignment_matrix
end

function longterm_cluster_IA_rates(channel, network, partition)
    I = get_no_BSs(network); K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    assignment = get_assignment(network)

    rates = zeros(Float64, K)

    for block in partition.blocks
        intercluster_interferers = setdiff(1:I, block.elements)
        for i in block.elements; for k in served_MS_ids(i, assignment)
            desired_power = channel.large_scale_fading_factor[k,i]^2*Ps[i]
            int_noise_power = sigma2s[k]
            for j in intercluster_interferers
                int_noise_power += channel.large_scale_fading_factor[k,j]^2*Ps[j]
            end
            rho = desired_power/int_noise_power

            rates[k] = 0.5log(1 + 2rho)
        end; end
    end

    return rates
end

include("BnBClustering.jl")
include("Chen2014_LinearObjClustering.jl")
include("ExhaustiveSearchClustering.jl")
include("GrandCoalitionClustering.jl")
include("NoClustering.jl")
include("RandomClustering.jl")
