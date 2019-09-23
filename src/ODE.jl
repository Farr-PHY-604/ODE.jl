module ODE

"""    stepsize_adjust(y, yerr, h, order; rtol=1e-8, atol=1e-8)

Returns `(hnew, retry_step)` giving a stepsize that will try to achieve the
given absolute and relative tolerances on the next step (or a retry of this
one).

`retry_step` is `true` if the error exceeds the bound by more than a factor of
two, and `false` otherwise.
"""
function stepsize_adjust(y, yerr, h, order; rtol=1e-8, atol=1e-8)
    scale_max = zero(y[1])
    for i in 1:length(y)
        scale = yerr[i]/(atol + rtol*abs(y[i]))
        scale_max = max(scale_max, scale)
    end

    retry_step = (scale_max > 2 ? true : false)

    hnew = h/(scale_max)^(1/(order+1))

    if hnew > 2*h
        hnew = 2*h
    elseif hnew < 1/2*h
        hnew = 1/2*h
    end

    (hnew, retry_step)
end

end # module
