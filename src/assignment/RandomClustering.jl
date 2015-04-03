# A random set partition is found and returned as the cluster. If the given
# partition has utility == -Inf, we keep searching. This can happen when IA
# infeasible clusters have -Inf utility.
function RandomClustering(channel, network)
    I = get_no_BSs(network); K = get_no_MSs(network)
    ds = get_no_streams(network); max_d = maximum(ds)

    aux_params = get_aux_assignment_params(network)
    @defaultize_param! aux_params "RandomClustering:max_iters" 1_000

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)

    random_partition = Partition()
    random_a = Array(Int, I)
    utilities = zeros(Float64, K, max_d)
    alphas = Array(Float64, K)

    # Find a set partition whose utility is not -Inf
    no_iters = 0; no_utility_calculations = 0
    while no_iters <= aux_params["RandomClustering:max_iters"]
        no_iters += 1

        # Get random partition by finding random rgs
        random_a = random_restricted_growth_string(I)
        random_partition = Partition(random_a)
        utilities, alphas, _ = longterm_utilities(channel, network, Partition(random_a))
        no_utility_calculations += 1

        if sum(utilities) > -Inf
            break
        end
    end
    if no_iters == aux_params["RandomClustering:max_iters"]
        Lumberjack.warn("Max iterations reached for RandomClustering. This probably means that an IA infeasible coalition structure was chosen.")
    end

    Lumberjack.info("RandomClustering finished.",
        { :sum_utility => sum(utilities),
          :a => random_a,
          :alphas => alphas  }
    )

    # Store alphas as user priorities for precoding, if desired
    if aux_params["apply_overhead_prelog"]
        set_user_priorities!(network, alphas)
    end

    # Store cluster assignment together with existing cell assignment
    temp_cell_assignment = get_assignment(network)
    network.assignment = Assignment(temp_cell_assignment.cell_assignment, cluster_assignment_matrix(network, random_partition))

    # Return results
    results = AssignmentResults()
    results["utilities"] = utilities
    results["a"] = random_a
    results["alphas"] = alphas
    results["no_iters"] = no_iters
    results["no_utility_calculations"] = no_utility_calculations
    return results
end

# Generate random restricted growth string. See Algorithm H in TAoCP 7.2.1.5.
function random_restricted_growth_string(I)
    a = Array(Int, I); b = Array(Int, I)
    for i = 1:I
        if i == 1
            b[1] = 1
            a[1] = 0
        else
            b[i] = 1 + maximum(a[1:(i-1)])
            a[i] = ifloor((b[i] + 1)*rand())
        end
    end
    return a
end
