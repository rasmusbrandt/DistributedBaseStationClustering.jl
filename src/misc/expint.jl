const expint_poly = Polynomials.Poly([
     8.267661952366478e+00,
    -7.773807325735529e-01,
    -3.012432892762715e-01,
    -7.811863559248197e-02,
    -1.019573529845792e-02,
    -6.973790859534190e-04,
    -2.569498322115933e-05,
    -4.819538452140960e-07,
    -3.602693626336023e-09
])

function expint(x::Float64)
    y = 0.
    eps1 = eps()

    if Polynomials.polyval(expint_poly, x) >= 0.
        # series expansion
        y = -0.57721566490153286061 - log(x)
        j = 1
        pterm = x
        term = x

        while abs(term) > eps1
            y = y + term
            j += 1
            pterm = -x*pterm/j
            term = pterm/j
        end
    else
        # continued fraction
        n = 1
        am2 = 0.
        bm2 = 1.
        am1 = 1.
        bm1 = x
        f = am1/bm1
        oldf = Inf
        j = 2

        while abs(f - oldf) > 100*eps1*abs(f)
            alpha = n - 1 + j/2

            a = am1 + alpha*am2
            b = bm1 + alpha*bm2

            am2 = am1/b
            bm2 = bm1/b
            am1 = a/b
            bm1 = 1.
       
            f = am1
            j += 1

            alpha = (j - 1)/2
            beta = x
            a = beta*am1 + alpha*am2
            b = beta*bm1 + alpha*bm2
            am2 = am1/b
            bm2 = bm1/b
            am1 = a/b
            bm1 = 1.
            oldf = f
            f = am1
            j += 1
        end

        y = exp(-x)*f
    end

    return y
end
