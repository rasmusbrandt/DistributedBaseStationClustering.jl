##########################################################################
# Random base station clustering.
#
# A random set partition is found and returned as the cluster. If the given
# partition has utility == -Inf, we keep searching. This can happen when IA
# infeasible clusters have -Inf utility.
function RandomClustering(channel, network)
    I = get_num_BSs(network); K = get_num_MSs(network)
    ds = get_num_streams(network); max_d = maximum(ds)

    aux_params = get_aux_assignment_params(network)
    @defaultize_param! aux_params "RandomClustering:max_iters" 1_000

    # Perform cell selection
    LargeScaleFadingCellAssignment!(channel, network)

    # Local RNG
    local_rng = MersenneTwister(round(Int, time()))

    random_partition = Partition(collect(0:(I-1))) # start with non-cooperative, if the loop below fails to find an IA feasible solution
    random_a = Array(Int, I)
    utilities = zeros(Float64, K, max_d)
    alphas = Array(Float64, K)

    # Find a set partition whose utility is not -Inf
    num_iters = 0
    num_utility_calculations = 0
    num_longterm_rate_calculations = 0
    while num_iters <= aux_params["RandomClustering:max_iters"]
        num_iters += 1

        # Get random partition by finding random rgs
        random_a = random_restricted_growth_string(I, local_rng)
        random_partition = Partition(random_a)
        utilities, alphas, _ = longterm_utilities(channel, network, Partition(random_a))
        num_utility_calculations += K
        num_longterm_rate_calculations += K

        if sum(utilities) > -Inf
            break
        end
    end
    if num_iters == aux_params["RandomClustering:max_iters"]
        Lumberjack.warn("Max iterations reached for RandomClustering. This probably means that an IA infeasible coalition structure was chosen.")
    end

    Lumberjack.info("RandomClustering finished.",
        @compat Dict(
            :sum_utility => sum(utilities),
            :a => random_a))

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
    results["num_clusters"] = 1 + maximum(random_a)
    results["num_iters"] = num_iters
    results["num_utility_calculations"] = num_utility_calculations
    results["num_longterm_rate_calculations"] = num_longterm_rate_calculations
    return results
end

# Generate random restricted growth string. See Algorithm H in TAoCP 7.2.1.5.
function random_restricted_growth_string(I, local_rng)
    a = Array(Int, I); b = Array(Int, I)
    for i = 1:I
        if i == 1
            b[1] = 1
            a[1] = 0
        else
            b[i] = 1 + maximum(a[1:(i-1)])
            a[i] = floor(Int, (b[i] + 1)*rand(local_rng))
        end
    end
    return a
end
