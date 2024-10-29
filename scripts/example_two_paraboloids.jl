# # Two Paraboloids with Constraints
# Import global pre-amble:
include("preamble.jl")
# Import callback type:
include("info_gathering_callback.jl")
#=
The following constrained problem is due to Gebken et. al.[^1]
\cite{gebkenDescentMethodEquality2018}.
\begin{equation}
  \begin{aligned}
    &\min_{\symbf{x}\in \mathbb{R}^2}
      \begin{bmatrix}
        (x_1 - 2)^2 + (x_2 - 1)^2 \\
        (x_1 - 2)^2 + (x_2 + 1)^2
      \end{bmatrix}
    \quad\text{s.t.}\quad
    g(\symbf{x}) = 1 - x_1^2 - x_2^2 \le 0.
  \end{aligned}
  \label{eqn:bennets_problem}
  \tag{Ex. 1}
\end{equation}
The feasible set of \eqref{eqn:bennets_problem} is $\mathbb{R}^2$ without the
interior of the unit ball.
The Pareto critical set is the line connecting $[2, -1]$ and $[2, 1]$ 
and the left boundary of the unit ball:
\begin{equation*}
  \mathcal P_c =
  \left\{
    \begin{bmatrix}
      2\\
      s
    \end{bmatrix}:
    s\in[-1,1]
  \right\}
  \bigcup
  \left\{
    \begin{bmatrix}
      \cos t\\
      \sin t
    \end{bmatrix}
    :
    t \in [\pi - \theta, \pi + \theta], \theta = \arctan\left(\frac{1}{2}\right)
  \right\}.
\end{equation*}
=#
#%%#src
# ## Problem Functions
# Before defining the actual objective and constraint functions, we set up
# helper types that allow for “simpler” counting of evaluation calls:

struct CountedFunction{F} <: Function
    wrapped_func :: F
    counter :: Base.RefValue{Int}
end
"""
    CountedFunction(func)

Wrap `func` and return a new `CountedFunction`.
"""
function CountedFunction(to_be_wrapped::F) where F
    return CountedFunction(to_be_wrapped, Ref(0))
end

# Calling a `CountedFunction` increases the internal counter and forwards
# to the wrapped function.
function (cfunc::CountedFunction)(args...; kwargs...)
    cfunc.counter[] += 1
    return cfunc.wrapped_func(args...; kwargs...)
end
# Finally, a helper function to reset the internal counter:
reset!(cfunc::CountedFunction)=(cfunc.counter[] = 0)
#%%
# Now for the actual objective functions.

# We define them as scalar valued functions first, so as to easier get Hessians
# later:
function objective_func1(x)
    return (x[1] - 2)^2 + (x[2] - 1)^2
end
function objective_func2(x)
    return (x[1] - 2)^2 + (x[2] + 1)^2
end
# Make them “counted”:
counted_objective_func1 = CountedFunction(objective_func1)
counted_objective_func2 = CountedFunction(objective_func2)

# The optimizer (or rather `MutableMOP`) wants vector-valued functions.
# For the objectives, we provide `counted_objectives`:
counted_objectives = function (x)
    return [
        counted_objective_func1(x);
        counted_objective_func2(x);
    ]
end

# In case of Taylor polynomial models, the gradient matrix is taken
# as the transposed, finite-difference Jacobian.
counted_objective_grads = CountedFunction(
    x -> transpose(FD.finite_difference_jacobian(counted_objectives, x))
)
# Likewise, the Hessian tensor is approximated using finite-differences:
counted_objective_hessians = CountedFunction(
    x -> begin
        H = zeros(2, 2, 2)
        H[:, :, 1] .= FD.finite_difference_hessian(counted_objective_func1, x)
        H[:, :, 2] .= FD.finite_difference_hessian(counted_objective_func2, x)
        return H
    end
)
# Setting up the (one-dimensional) constraint function works similarly:
function vec_constr_func(x)
    return [ 1 - sum( x.^2 ), ]
end
counted_constrs = CountedFunction(vec_constr_func)
## Gradient:
counted_constr_grads = CountedFunction(
    x -> transpose(FD.finite_difference_jacobian(counted_constrs, x))
)
## Hessian:
counted_constr_hessians = CountedFunction(
    x -> FD.finite_difference_hessian(_x -> counted_constrs(_x)[1], x)
)

# Let's also define a helper function to reset all counters:
function reset_all_counters!()
    funcs = [
        counted_objective_func1, counted_objective_func2, counted_objective_grads,
        counted_objective_hessians,
        counted_constrs, counted_constr_grads, counted_constr_hessians
    ]
    for counted_func in funcs
        reset!(counted_func)
    end
    return nothing
end
# The final optimization problem is built at a later point using the following helper:
function prepare_mop(model_cfg, constr_type=:ineq)
    mop = C.MutableMOP(; num_vars = 2)
    ## add the vector-valued objective
    C.add_objectives!(
        mop, counted_objectives, model_cfg;
        dim_out=2, func_iip=false,  # `func_iip=false` because function is out-of-place, not in-place
        grads = counted_objective_grads, grads_iip = false,
        hessians = counted_objective_hessians, hessians_iip = false
    )
    ## add constraints
    ## … inequality or equality constraints?
    add_constrs! = constr_type == :ineq ? C.add_nl_ineq_constraints! : C.add_nl_eq_constraints!
    ## … add functions to `mop`
    add_constrs!(
        mop, counted_constrs, model_cfg;
        dim_out=1, func_iip=false,
        grads = counted_constr_grads, grads_iip = false,
        hessians = counted_constr_hessians, hessians_iip = false
    )
    return mop
end
# Which of the above functions are actually required depends on the surragate model configuration.

# ## Experiments
# ### Comparing Function Calls

# Options for the trust-region algorithm:
algo_opts = C.AlgorithmOptions(;
    log_level = Logging.Debug,
    stop_delta_min = 1e-9,
    nu_accept = 0.01,
    nu_success = 0.9,
    kappa_theta = 1e-4,
    gamma_shrink_much = 0.25,
    gamma_shrink = .5,
    gamma_grow = 2.0,
    eps_crit = 1e-4,
    eps_theta = 1e-6,
    psi_theta = 2.0,
    crit_B = 1000,
    crit_M = 3000,
    c_delta = 0.99,
    c_mu = 100,
    mu = 0.01,
    stop_xtol_rel = eps(Float64) * 100,
    trial_mode = Val(:max_diff),
    #trial_mode = Val(:max_rho),
    trial_update = Val(:stepsize),
    step_config = C.SteepestDescentConfig(;
        backtracking_mode = Val(:max),
        rhs_factor = 1e-4,
        normal_step_norm = 1,
        descent_step_norm = Inf,
    )
)
#%%#src

# Next, define the different surrogate models that we want to compare.
# Store surrogate configuration and other settings in `experiment_configs`.

## RBF settings:
kernel = C.CubicKernel()
max_points = 6

## algorithmic parameters to be varied
nu_success_a = .45
nu_success_b = .9

## complete settings dict:
experiment_configs = OrderedDict(
    "rbf_a" => (; x0 = [-2.0, 0.5], nu_success = nu_success_a, model_cfg = C.RBFConfig(; kernel, max_points)),
    "tp1_a" => (; x0 = [-2.0, 0.5], nu_success = nu_success_a, model_cfg = C.TaylorPolynomialConfig(; degree=1)),
    "tp2_a" => (; x0 = [-2.0, 0.5], nu_success = nu_success_a, model_cfg = C.TaylorPolynomialConfig(; degree=2)),
    "rbf_b" => (; x0 = [-2.0, -0.5], nu_success = nu_success_b, model_cfg = C.RBFConfig(; kernel, max_points)),
    "tp1_b" => (; x0 = [-2.0, -0.5], nu_success = nu_success_b, model_cfg = C.TaylorPolynomialConfig(; degree=1)),
    "tp2_b" => (; x0 = [-2.0, -0.5], nu_success = nu_success_b, model_cfg = C.TaylorPolynomialConfig(; degree=2)),
    "rbf_aa" => (; x0 = [-2.0, 0.0], nu_success = nu_success_a, model_cfg = C.RBFConfig(; kernel, max_points)),
)
# The results are stored in a dict-of-dicts, and the keys stem from
# `experiment_configs`.
results = OrderedDict{String, Any}()

# We iterate the keys and surrogate configurations to set up the optimization problems,
# run them and store some performance indicators:
for (model_key, cfg) in pairs(experiment_configs)
    global algo_opts
    @unpack model_cfg, x0, nu_success = cfg
    ## reset all function counters
    reset_all_counters!()

    @reset algo_opts.nu_success = nu_success
    
    ## initialize a 2D optimization problem
    mop = prepare_mop(model_cfg, :ineq)
   
    ## initialize empty callback
    user_callback = InfoGatheringCallback()
    ## run optimization
    ret_obj = C.optimize(mop, x0; algo_opts, user_callback)
    ## process and store results for plotting
    results[model_key] = Dict(
        "x0" => copy(x0),
        "xfin" => C.opt_vars(ret_obj),
        "x" => reduce(hcat, user_callback.x),
        "nc_primal" => max(
            counted_objective_func1.counter[],
            counted_objective_func2.counter[],
            counted_constrs.counter[]
        ),
        "nc_diff" => max(
            counted_objective_grads.counter[],
            counted_constr_grads.counter[],
        ),
        "nc_hess" => max(
            counted_objective_hessians.counter[],
            counted_constr_hessians.counter[],
        ),
        "num_its" => isempty(user_callback.x) ? 0 : length(user_callback.x) - 1
    )
end
#%%#src
# #### Plot Figure 1
# First, enable our custom theme:
set_theme!(CUSTOM_THEME)

# Initialize the figure:
fig = Figure(; size=(483, 310), figure_padding=2f0)

# Add title labels:
Label(
    fig[1, 1:2],
    rich(
        "Trajectories with thresholds ", 
        rich("ν", subscript("(1)"), " = $(nu_success_a)"; color=:black),
        " and ",
        rich("ν", subscript("(1)"), " = $(nu_success_b)"; color=:red),
        "."
        ;
        font = :bold
    );
    fontsize = 12f0
)
Label(
    fig[2, 1:2],
    "(№ iterations, № function calls, № gradient calls, № Hessian calls)";
)

ax = Axis(
    fig[3, 1];
    ##aspect = DataAspect(),
    autolimitaspect = 1,
    xlabel = L"x_1",
    ylabel = L"x_2",
)

# Plot infeasible set:
poly!(
    ax, [Point2f(cos(phi), sin(phi)) for phi in range(0, 2*π, 100)];
    color = (:yellow, .5f0),
    label = "infeasible"
)

# Plot Pareto-Set:
## helper
function plot_unit_arcs!(ax, angs...; 
    color=:black, linewidth=3f0, single_label = nothing, kwargs...
)
    for (i,ang) in enumerate(angs)
        label = i == 1 ? single_label : nothing
        R = range(ang - atan(0.5), ang + atan(0.5), 100)
        lines!(
            ax, cos.(R), sin.(R); 
            color, linewidth, label, kwargs...
        )
    end
end
## plot the Pareto arc
plot_unit_arcs!(ax, π; label="Pareto Set")
## and the Pareto line
lines!(
    ax, [(2.0, -1.0), (2.0, 1.0)]; 
    color = :black,
    linewidth = 3f0,
)

# For the trajectories, plotting properties are defined in a Dict:
plot_props = Dict(
    "rbf_a" => (label_prefix = "RBF", color = Makie.wong_colors()[1], marker=:circle),
    "rbf_b" => (label_prefix = "RBF", color = Makie.wong_colors()[1], marker=:circle),
    "tp1_a" => (label_prefix = "TP1", color = Makie.wong_colors()[2], marker=:rect),
    "tp1_b" => (label_prefix = "TP1", color = Makie.wong_colors()[2], marker=:rect),
    "tp2_a" => (label_prefix = "TP2", color = Makie.wong_colors()[3], marker=:diamond),
    "tp2_b" => (label_prefix = "TP2", color = Makie.wong_colors()[3], marker=:diamond),
    "rbf_aa" => (label_prefix = "RBF", color = Makie.wong_colors()[4], marker=:circle),
)

# For padding in mono-spaced labels, compute the number of digits for all 
# relevant data points:
_result_numbers = res_dict -> (
    res_dict["num_its"], res_dict["nc_primal"], res_dict["nc_diff"], res_dict["nc_hess"])
label_padding = zeros(Int, 4)
for res_dict in values(results)
    global label_padding
    label_padding = max.(label_padding, ndigits.(_result_numbers(res_dict)))
end 

# Loop through results and make trajetory plots:
for (model_key, res_dict) in pairs(results)
    global label_padding
    model_props = plot_props[model_key]
    
    @unpack label_prefix, color, marker = model_props
    num_its, nc_primal, nc_diff, nc_hess = _result_numbers(res_dict)
    label = label_prefix * 
        "($(num_its), $(nc_primal), $(nc_diff), $(nc_hess))"
    
    label = rich(
        label_prefix, 
        " ($(lpad("$(num_its)", label_padding[1], " "))",
        ",$(lpad("$(nc_primal)", label_padding[2], " "))",
        ",$(lpad("$(nc_diff)", label_padding[3], " "))",
        ",$(lpad("$(nc_hess)", label_padding[4], " ")))",
        font="Latin Modern Mono Light"
    )

    strokecolor = endswith(model_key, "a") ? :black : :red
    scatterlines!(
        ax, res_dict["x"]; 
        label, color, marker, strokecolor,
    )
    scatter!(ax, Point2f(res_dict["xfin"]...);
        color, marker, strokecolor,
        markersize = 8.5f0,
        strokewidth = 1.5f0,
    )
end

# Make legend:
Legend(
    fig[3, 2], ax; 
    framevisible = true, 
)
# Tinker with layout:
rowgap!(fig.layout, 10f0)
colgap!(fig.layout, 10f0)

# Save the plot:
save(
    ensure_path(joinpath(PLOTS_PATH, "two_paraboloids_ineq.png")), fig; 
    pt_per_unit=1, px_per_unit=5
)
fig
#%%#src
# ### Equality Constraint
# If we use the constraint function as an equality constraint instead of 
# an inequality constraint, the left-most part of the unit sphere is Pareto-critical,
# as is the right-most part.
# But the connecting line of individual minima is not critical anymore.
#
# We test this using RBF surrogates.
# First, reset some of the global options:
@reset algo_opts.nu_success = .9
@reset algo_opts.max_iter = 1000
@reset algo_opts.log_level = Logging.Debug

# Then set-up the MOP:
mop = prepare_mop(C.RBFConfig(; kernel=C.CubicKernel(), max_points=6), :eq)
# We compare different initial trust-region sizes, starting points and 
# compatibility paramaters:
configs_cubics = (
    (; x0 = 1.5 .* [cos(8/9 * π), sin(8/9 * π)], delta_init=1.0, c_delta=.99, ),
    (; x0 = 1.25 .* [cos(8.5/9 * π), sin(8.5/9 * π)], delta_init=1.0, c_delta=.5, ),
    (; x0 = 1.25 .* [cos(9.5/9 * π), sin(9.5/9 * π)], delta_init=0.1, c_delta=.5, ),
    (; x0 = 1.5 .* [cos(10/9 * π), sin(10/9 * π)], delta_init=1.0, c_delta=.5, ),
)
# Run all the experiments:
results_cubics = Any[]
for cfg in configs_cubics
    global algo_opts
    ## set algo_opts
    @reset algo_opts.delta_init = cfg.delta_init
    @reset algo_opts.c_delta = cfg.c_delta
    ## ensure correct call numbers
    reset_all_counters!()
    ## make a new callback
    user_callback = InfoGatheringCallback()
    ## run
    ret_cubic = C.optimize(
        mop, cfg.x0;
        algo_opts, user_callback
    )
    ## Extract number of iterations and function calls:
    num_its = length(user_callback.x)
    nc_objf = max(counted_objective_func1.counter[], counted_objective_func2.counter[])
    nc_constr = counted_constrs.counter[]
    push!(
        results_cubics,
        (; 
            cback = deepcopy(user_callback), ret = deepcopy(ret_cubic), 
            cfg, num_its, nc_objf, nc_constr
        )
    )
end

#%%#src
# #### Plot Figure 2
# Set up figure …
fig2 = Figure(; size = (483, 320))
# … and sub-layouts for tricky placing:
top_layout = GridLayout()
bottom_layout = GridLayout()

# Axis and title:
ax1 = bottom_layout[1, 1] = Axis(
    fig2;
    limits=((-1.75, 2.75), nothing),
    xticks = -2:1:3,
    yticks = -1:1:1,
    autolimitaspect=1f0,
    xlabel = L"x_1",
    ylabel = L"x_2",
)

top_layout[1, 1:2] = Label(
    fig2,
    rich("Trajectories with Cubic Kernel."; font=:bold);
    fontsize=12f0
)

# Plot Pareto-Set:
## unit sphere for visualization
unit_range = range(0, 2*π, 100)
lines!(
    cos.(unit_range), sin.(unit_range); 
    color=:gray70, linestyle=:dash, 
    label = L"\mathcal S^1",
)
plot_unit_arcs!(ax1, 0, π; single_label="Pareto Set")

# Plot iterations
function plot_rbf_trajectory!(ax, cback; kwargs...)
    x = reduce(hcat, cback.x)
    ## make markers dependent on iteration type
    marker = Symbol[]
    markersize = Float32[]
    for it_stat in cback.it_status
        it_class = it_stat.iteration_classification
        m, ms = if it_class == C.IT_RESTORATION 
            :star6, 9f0
        else
            if it_class == C.IT_THETA_STEP 
                :rect, 7f0
            else
                :circle, 7f0
            end
        end
        push!(marker, m)
        push!(markersize, ms)
    end
    l = lines!(ax, x; kwargs...)
    s = scatter!(ax, x; marker, markersize)
    return l, s
end
for (i, res) in enumerate(results_cubics)
    label = rich(
        rich("Δ", subscript("(0)"), " = $(res.cfg.delta_init), "),
        rich(
            rich("c"; font="Latin Modern Mono Light"), 
                subscript("Δ"), " = $(res.cfg.c_delta)")
    )
    color = Makie.wong_colors()[i]
    plot_rbf_trajectory!(ax1, res.cback; label, color)
    xfin = Point2(C.opt_vars(res.ret)...)
    text!(
        ax1, xfin; 
        text = "($(res.num_its), $(res.nc_objf), $(res.nc_constr))",
        color,
        align = xfin[1] < 0 ? (:left, :center) : (:right, :center),
        offset = xfin[1] < 0 ? (8f0, 0f0) : (-8f0, 0f0),
        fontsize = 12f0
    )
end
# Legend for trajectories:
axislegend(; 
    ax1, position = :rb
)

# Manual Legend
top_layout[2,1] = Legend(
    fig2,
    [
        MarkerElement(color=:black, marker=:rect),
        MarkerElement(color=:black, marker=:star6),
        MarkerElement(color=:black, marker=:circle),
    ],
    [
        L"$θ$-Step",
        "Restoration",
        "Other"
    ];
    orientation = :horizontal,
    padding = (3f0, 3f0, 3f0, 3f0),
    patchsize = (5f0, 10f0),
    halign = :left
)
top_layout[2,2] = Label(
    fig2,
    "(№ iterations, № objective calls, № constraint calls)";
    halign = :left
)

# Place sub-layouts:
fig2.layout[1, 1] = top_layout
fig2.layout[2, 1] = bottom_layout
# Tweak sizes:
colsize!(bottom_layout, 1, Fixed(428))
rowgap!(top_layout, 5f0)
rowgap!(fig2.layout, 10f0)
colgap!(top_layout, 0)
# Save and show:
save(
    ensure_path(joinpath(PLOTS_PATH, "two_paraboloids_eq.png")), fig2;
    pt_per_unit=1,
    px_per_unit=5
)

fig2
#%%#src
#=
[^1]: B. Gebken, S. Peitz, and M. Dellnitz,
    “A Descent Method for Equality and Inequality Constrained 
    Multiobjective Optimization Problems,”
    Numerical and Evolutionary Optimization – 
    NEO 2017, 2018, 
    Available: https://ris.uni-paderborn.de/record/8750
=#