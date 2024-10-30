# # MW3 Example

## Imports
# General preamble for examples:
include("preamble.jl")

# A callback to gather information from Compromise:
include("info_gathering_callback.jl")

import NLopt
import HaltonSequences: HaltonPoint

# ## MW3 Functions
const sqrt2 = sqrt(2)
d(x) = 1 + 2 * (x[2] + (x[1] - 0.5)^2 - 1)^2 + 
    2 * (x[3] + (x[2] - 0.5)^2 - 1)^2

f1(x) = x[1]
f2(x) = d(x) * (1 - (f1(x) / d(x)))
f(x) = [f1(x); f2(x)]

l(f1x, f2x) = sqrt2 * (f2x - f1x)
c1(f1x, f2x) = +f1x + f2x - 1.05 - 0.45 * (sin(3/4 * π * l(f1x, f2x)))^6
c2(f1x, f2x) = -f1x - f2x + 0.85 + 0.30 * (sin(3/4 * π * l(f1x, f2x)))^2
c1(x) = c1(f1(x), f2(x))
c2(x) = c2(f1(x), f2(x))
c(x) = begin
    f1x = f1(x)
    f2x = f2(x)
    return [
        c1(f1x, f2x),
        c2(f1x, f2x),
    ]
end

# Utility functions for plotting constraints in objective space:
cmax(f1x, f2x) = max(c1(f1x, f2x), c2(f1x, f2x))
copt(f1x, f2x) = max(c2(f1x, f2x), 1 - f2x - f1x)
cfeas(f1x, f2x) = max(copt(f1x, f2x), c1(f1x, f2x))

# Lower and upper variable bounds:
lb = zeros(3)
ub = ones(3)

max_func_calls = 100
#%%
# Find starting points
num_points = 15
X0vec = Vector{Float64}[]
for x0 in HaltonPoint(3)
    if f2(x0) < 1.5
        push!(X0vec, x0)
        if length(X0vec) >= num_points
            break
        end
    end
end
X0 = reduce(hcat, X0vec)
F0 = mapreduce(f, hcat, X0vec)
#%%
# ## Compromise
function compromise_trajectory(x0; kwargs...)
    global lb, ub, f, c
    global max_func_calls
    
    rbf_config = C.RBFConfig()
    algo_opts = C.AlgorithmOptions(; 
        stop_delta_min = 1e-10,
        log_level = Debug,
        backtrack_in_crit_routine=true,
        kwargs...
    )
    
    mop = C.MutableMOP(; num_vars=3, lb, ub)
    C.add_objectives!(mop, f, rbf_config; 
        dim_out = 2, func_iip=false, max_func_calls)
    C.add_nl_ineq_constraints!(mop, c, rbf_config; 
        dim_out = 2, func_iip=false, max_func_calls)

    user_callback = InfoGatheringCallback()
    ret = C.optimize(mop, x0; user_callback, algo_opts)
    
    X = hcat(reduce(hcat, user_callback.x), C.opt_vars(ret))
    F = mapreduce(f, hcat, eachcol(X))

    return X, F
end

function compromise_trajectories(X0; kwargs...)
    X = Vector{Matrix{Float64}}()
    F = Vector{Matrix{Float64}}()
    for x0 = eachcol(X0)
        _X, _F = compromise_trajectory(x0; kwargs...)
        push!(X, _X)
        push!(F, _F)
    end
    return X, F
end

# ## Weighted Sum (NLopt)
function nlopt_trajectory(x0; 
    ctol=1e-10,
)
    global f, c, lb, ub
    global max_func_calls
    opt = NLopt.Opt(:LN_COBYLA, 3)
    NLopt.lower_bounds!(opt, lb)
    NLopt.upper_bounds!(opt, ub)

    nlopt_objective = function(x, dfx)
        @assert isempty(dfx) "Derivative-free algorithms only."
        return f1(x) + f2(x)
    end
    NLopt.min_objective!(opt, nlopt_objective)

    nlopt_constraint = function(cx, x, Dcx)
        @assert isempty(Dcx) "Derivative-free algorithms only."
        cx .= c(x)
        return cx
    end
    NLopt.inequality_constraint!(opt, nlopt_constraint, fill(ctol, 2))

    NLopt.maxeval!(opt, max_func_calls)
    NLopt.xtol_abs!(opt, 1e-10)
    NLopt.xtol_rel!(opt, 0)
    NLopt.ftol_abs!(opt, 0)
    NLopt.ftol_rel!(opt, 0)
    (fopt, xopt, ret) = NLopt.optimize(opt, x0)

    f0 = f(x0)
    fopt = f(xopt)
    return hcat(x0, xopt), hcat(f0, fopt)
end
function nlopt_trajectories(X0; kwargs...)
    X = Vector{Matrix{Float64}}()
    F = Vector{Matrix{Float64}}()
    for x0 = eachcol(X0)
        _X, _F = nlopt_trajectory(x0)
        push!(X, _X)
        push!(F, _F)
    end
    return X, F
end
#%%
# ## Run the Experiments
Xn, Fn = nlopt_trajectories(X0);
Xc, Fc = compromise_trajectories(X0;
    trial_mode = Val(:max_diff),
    step_config = C.SteepestDescentConfig(;
        backtracking_mode = Val(:max)
    )
);
Xc_, Fc_ = compromise_trajectories(X0;
    trial_mode = Val(:min_rho),
    step_config = C.SteepestDescentConfig(;
        backtracking_mode = Val(:all),
        normalize_gradients = true,
    ),
);
#%%
# ## Plotting
set_theme!(CUSTOM_THEME; Scatter = (; markersize=10f0))

F1 = LinRange(0.0, 1.0, 200)
F2 = copy(F1)

fig = Figure(; size = (480, 250), figure_padding=2f0)
ax = Axis(
    fig[1,1];
    title="Multi-Objective Optimization vs Weighted Sum",
    xlabel=L"f_1",
    ylabel=L"f_2",
    xlabelpadding=2f0,
    ylabelpadding=2f0
)
# Plot all image space vectors satisfying constraints:
contourf!(
    ax, F1, F2, cmax; 
    levels = [0,], extendlow=(:black, .2)
)
# Plot only those image space vectors that have a feasible pre-image:
contourf!(
    ax, F1, F2, cfeas; 
    levels = [0,], extendlow=(:black, .2)
)
# Plot the Pareto front:
contour!(
    ax, F1, F2, copt; 
    levels = [0,], color=:black,
    linewidth = 3f0,
    label="PS"
)

# Starting points:
for (i,y) in enumerate(eachcol(F0))
    scatter!(ax, (y...); color=Cycled(i), marker=:circle)
end
# Helper for entire trajectories:
function plot_trajectories!(
    ax, F;
    marker=:circle,
    linestyle=:solid
)
    for (i,Y) in enumerate(F)
        lines!(ax, Y; alpha=.35, color=Cycled(i), linestyle)
        _Y = Y[:, [end,]]
        scatter!(ax, _Y; color=Cycled(i), marker)
    end
end
# Plot trajectories using different symbols
plot_trajectories!(ax, Fc; marker = :star5)
plot_trajectories!(ax, Fc_; marker = :rect)
plot_trajectories!(ax, Fn; marker = :cross, linestyle=:dash)

# Hack a legend
font = "Latin Modern Mono Light"
lines!(ax, zeros(2,0); color=(:black, .4), label="feasible", linewidth=10f0)
scatter!(ax, zeros(2,0);
    marker=:star5, color=:white, label=rich("Compromise std."; font))
scatter!(ax, zeros(2,0);
    marker=:rect, color=:white, label=rich("Compromise mod."; font))
scatter!(ax, zeros(2,0); 
    marker=:cross, color=:white, label=rich("WS COBYLA"; font))

axislegend(ax; position=(.1, .1))
save(
    ensure_path(joinpath(PLOTS_PATH, "mw3_trajectories.pdf")), fig; 
    pt_per_unit=1, px_per_unit=5
)

fig