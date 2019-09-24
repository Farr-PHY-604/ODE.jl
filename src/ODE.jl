module ODE

export evolve_ode, Midpoint

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

struct Stepper
    order::Int
    step
end

"""    evolve_ode(stepper, f_rhs, x0, y0, x1; h_init, atol, rtol)

Returns the state of the ODE represented by `f_rhs` evolved from its initial
condition `(x0, y0)` to `x1` using the single-step algorithm in `stepper`.

`h_init` is used to guess the first step.

The evolution is adjusted to control the the absolute and relative error in the
solution according to the given absolute and relative tolerances; steps may be
re-tried until an acceptable error is achieved.

The final step will land *exactly* on `x1`.
"""
function evolve_ode(stepper, f_rhs, x0, y0, x1; h_init=(x1-x0)/100, atol=1e-8, rtol=1e-8)
    x = x0
    y = y0
    h = h_init
    not_done = true

    while not_done
        if x + h > x1
            h = x1-x
            not_done = false
        end

        y_new, yerr = stepper.step(f_rhs, x, y, h)
        h_new, retry = stepsize_adjust(y_new, yerr, h, stepper.order; rtol=rtol, atol=atol)

        if retry
            h = h_new
            not_done = true
        else
            y = y_new
            x = x + h
            h = h_new
        end
    end

    y
end

function step_midpoint(f_rhs, x, y, h)
    dydx0 = f_rhs(x, y)
    h_half = h/2
    y_half = y .+ dydx0 .* h_half

    dydx_half = f_rhs(x+h_half, y_half)

    y_new = y .+ dydx_half .* h
    yerr = h .* abs.(dydx_half .- dydx0)

    (y_new, yerr)
end

Midpoint = Stepper(2, step_midpoint)


end # module
