# See description in IntraclusterWMMSEState for how to interpret
# the different rates.
function calculate_weighted_logdet_rates(state, alphas)
    K = length(state.V)
    ds = Int[ size(state.V[k], 2) for k = 1:K ]; max_d = maximum(ds)

    weighted_logdet_rates_full = zeros(Float64, K, max_d)
    weighted_logdet_rates_partial = zeros(Float64, K, max_d)
    weighted_logdet_rates_LB = zeros(Float64, K, max_d)
    for k = 1:K; if ds[k] > 0
        # W is p.d., so we should only get abs eigenvalues. Numerically we may
        # get some imaginary noise however. Also, numerically the eigenvalues
        # may be less than 1, so we need to handle that to not get negative
        # rates. (Same goes for E, but in reverse.)
        r_full = alphas[k]*log2(1./(min(1, abs(eigvals(state.E_full[k])))))
        r_partial = alphas[k]*log2(1./(min(1, abs(eigvals(state.E_partial[k])))))
        r_lower_bound = alphas[k]*log2(1./(min(1, abs(eigvals(state.E_LB[k])))))

        if ds[k] < max_d
            weighted_logdet_rates_full[k,:] = cat(1, r_full, zeros(Float64, max_d - ds[k]))
            weighted_logdet_rates_partial[k,:] = cat(1, r_partial, zeros(Float64, max_d - ds[k]))
            weighted_logdet_rates_LB[k,:] = cat(1, r_lower_bound, zeros(Float64, max_d - ds[k]))
        else
            weighted_logdet_rates_full[k,:] = r_full
            weighted_logdet_rates_partial[k,:] = r_partial
            weighted_logdet_rates_LB[k,:] = r_lower_bound
        end
    end; end

    return weighted_logdet_rates_full, weighted_logdet_rates_partial, weighted_logdet_rates_LB
end

##########################################################################
# Hermitian methods used in precoding methods
if VERSION >= v"0.4"
    macro hermdata(A)
        return :($A.data)
    end
else
    macro hermdata(A)
        return :($A.S)
    end
end

import Base: +, -, .*
+(A::Hermitian{Complex128}, B::Hermitian{Complex128}) = Hermitian(@hermdata(A) + @hermdata(B))
+(B::Matrix{Float64}, A::Hermitian{Complex128}) = +(A, B)
+(A::Hermitian{Complex128}, B::Matrix{Float64}) = @hermdata(A) + B
+(A::Hermitian{Complex128}, B::Matrix{Complex128}) = @hermdata(A) + B

-(A::Hermitian{Complex128}, B::Matrix{Complex128}) = @hermdata(A) - B
-(B::Array{Complex128, 2}, A::Hermitian{Complex128}) = -(A, B)
-(A::Hermitian{Complex128}, B::Matrix{Float64}) = @hermdata(A) - B

.*(a::Float64, B::Hermitian{Complex128}) = Hermitian(a.*@hermdata(B))

import Base.logdet, Base.diag
logdet(A::Hermitian{Complex128}) = logdet(@hermdata(A))
diag(A::Hermitian{Complex128}) = diag(@hermdata(A))
