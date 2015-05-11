expint(x::Real) = expint([x])[1]

function expint{T<:Real}(x::Vector{T})
    y = zeros(Float64, length(x))

    p = Polynomials.Poly([8.267661952366478e+00, -7.773807325735529e-01, -3.012432892762715e-01, -7.811863559248197e-02, -1.019573529845792e-02, -6.973790859534190e-04, -2.569498322115933e-05, -4.819538452140960e-07, -3.602693626336023e-09])
    polyv = Polynomials.polyval(p, x)

    # series expansion
    k = find(polyv .>= 0.)
    if !isempty(k)
        egamma = 0.57721566490153286061
        xk = x[k]
        yk = -egamma - log(xk)
        j = 1
        pterm = xk
        term = xk

        while any(abs(term) .> eps())
            yk = yk + term
            j = j + 1
            pterm = -xk.*pterm/j
            term = pterm/j
        end

       y[k] = yk
    end

    # continued fraction
    k = find(polyv .< 0.)
    if !isempty(k)
        n = 1
        xk = x[k]
        am2 = zeros(xk)
        bm2 = ones(xk)
        am1 = ones(xk)
        bm1 = xk;
        f = am1 ./ bm1
        oldf = fill(Inf,length(xk))
        j = 2

        while any(abs(f-oldf) .> (100*eps().*abs(f)))
            alpha = n-1+(j/2)

            a = am1 + alpha * am2
            b = bm1 + alpha * bm2

            am2 = am1 ./ b
            bm2 = bm1 ./ b
            am1 = a ./ b
            bm1 = 1
       
            f = am1
            j = j+1

            alpha = (j-1)/2
            beta = xk
            a = beta .* am1 + alpha * am2
            b = beta .* bm1 + alpha * bm2
            am2 = am1 ./ b
            bm2 = bm1 ./ b
            am1 = a ./ b
            bm1 = 1
            oldf = f
            f = am1
            j = j+1
        end

        y[k] = exp(-xk) .* f
    end

    return y
end
