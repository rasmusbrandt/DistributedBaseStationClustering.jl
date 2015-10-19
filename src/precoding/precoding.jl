# See description in IntraclusterWMMSEState for how to interpret
# the different rates.
function calculate_weighted_logdet_rates(state, prelogs_cluster_sdma, prelogs_network_sdma)
    K = length(state.E_cluster_sdma_full)
    ds = Int[ size(state.E_cluster_sdma_full[k], 2) for k = 1:K ]; max_d = maximum(ds)

    weighted_logdet_rates_cluster_sdma_full    = zeros(Float64, K, max_d)
    weighted_logdet_rates_network_sdma_full    = zeros(Float64, K, max_d)
    weighted_logdet_rates_network_sdma_partial = zeros(Float64, K, max_d)
    for k = 1:K; if ds[k] > 0
        weighted_logdet_rates_cluster_sdma_full[k,1:ds[k]]    = prelogs_cluster_sdma[k]*(-1)*log2(abs(eigvals(state.E_cluster_sdma_full[k])))
        weighted_logdet_rates_network_sdma_full[k,1:ds[k]]    = prelogs_network_sdma[k]*(-1)*log2(abs(eigvals(state.E_network_sdma_full[k])))
        weighted_logdet_rates_network_sdma_partial[k,1:ds[k]] = prelogs_network_sdma[k]*(-1)*log2(abs(eigvals(state.E_network_sdma_partial[k])))
    end; end

    return weighted_logdet_rates_cluster_sdma_full, weighted_logdet_rates_network_sdma_full, weighted_logdet_rates_network_sdma_partial
end

# This helper is needed in 0.4, since the diagonals must necessarily be reals.
function hermitianize!(M)
    n = Base.LinAlg.chksquare(M)
    for i = 1:n
        M[i,i] = real(M[i,i])
    end
    Hermitian(M)
end
