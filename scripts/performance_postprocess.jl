# # Performance Profiles

# In this script, we performe some analysis on the optimization results.
# To this end, we have to load the test problem CSV files into dataframes
## dependencies
include("preamble.jl")
## constrained problems
lib_constr = MOB.readlib(LIB_CONSTR_CSV_FPATH);
## unconstrained problems
lib_unconstr = MOB.readlib(LIB_UNCONSTR_CSV_FPATH);

# We assume that the CSV files have been created in such a way that they
# have columns for all the solvers as given by these arrays:
const SOLVERS_UNCONSTR = (
    COL_NAME_DFMO,
    COL_NAME_COMPROMISE_SEQ,
    COL_NAME_COMPROMISE_SET
)
const SOLVERS_CONSTR = (
    SOLVERS_UNCONSTR...,
    COL_NAME_COMPROMISE_PEN
)
# Most metrics need a combined Pareto-front.
# We compute row-wise combined fronts by doing non-dominance testing.
# Infeasible solutions are not considered.
viol_tol = 1e-3     # maximum constraint violation
MOB.combine_fronts!(
    lib_unconstr, SOLVERS_UNCONSTR;
    ## column to be added
    cname = "front_fpath",
    viol_tol,
    overwrite = false
)

MOB.combine_fronts!(
    lib_constr, SOLVERS_CONSTR;
    cname = "front_fpath",
    viol_tol, 
    overwrite = false
)

# Store the modified CSV files:
MOB.store!(lib_unconstr; outname=LIB_UNCONSTR_CSV_NAME)
MOB.store!(lib_constr; outname=LIB_CONSTR_CSV_NAME)

#%%
# ## Purity Performance
# The “purity” metric is best compared in pairs.
# This simple helper function generates solver column name pairs:

function pairwise_combinations(objects)
    T = eltype(objects)
    obj_pairs = Vector{Tuple{T, T}}()
    N = length(objects)
    for i = 1 : N
        obj_i = objects[i]
        for j = i+1 : N
            obj_j = objects[j]
            push!(obj_pairs, (obj_i, obj_j))
        end
    end
    return obj_pairs
end
# Generate the solver pair arrays:
PAIRS_UNCONSTR = pairwise_combinations(SOLVERS_UNCONSTR)
PAIRS_CONSTR = pairwise_combinations(SOLVERS_CONSTR)

# Before computing the performance and generating plots, I define a global property
# dict to have nicer labels etc.
solver_props = Dict(
    COL_NAME_DFMO => (; name="DFMO", color=Makie.wong_colors()[1], marker=:circle),
    COL_NAME_COMPROMISE_SEQ => (; name="CSEQ", color=Makie.wong_colors()[2], marker=:utriangle),
    COL_NAME_COMPROMISE_SET => (; name="CSET", color=Makie.wong_colors()[3], marker=:star5),
    COL_NAME_COMPROMISE_PEN => (; name="CPEN", color=Makie.wong_colors()[4], marker=:cross),
)

# The following function takes a tuple of solver column names (from `PAIRS_UNCONSTR`
# or `PAIRS_CONSTR`) and computes the purity metric.
# A new file is created based on the solver names and the purity results are stored.
function purity_for_solver_pair(
    results_df, solver_fpath_pair;
    ref_front_cname = "front_fpath",
    prefix = "", suffix = ""
)
    global solver_props, viol_tol
    ## unpack column names
    solverA_fpath, solverB_fpath = solver_fpath_pair
    ## get names
    solverA = solver_props[solverA_fpath].name
    solverB = solver_props[solverB_fpath].name
    ## make filename
    results_basename = prefix * "purity_" * solverA * "_" * solverB * suffix
    results_name = results_basename * ".jld2"
    ## compute performance
    @info "Computing Purity for $(solverA) & $(solverB)."
    purity_df = MOB.purity_with_reference(
        results_df, solver_fpath_pair, ref_front_cname;
        overwrite=false,
        results_name,
        viol_tol
    )
    DF.metadata!(purity_df, "results_basename", results_basename; style=:note)
    return purity_df
end

# When we have a purity dataframe, we still need to turn it into a performance profile.
# This is done with `MOB.performance_profile`.
# Afterwards, we can make the plots:

function plot_performance_profile(
    perf_df, solver_fpath_cnames = nothing; 
    title = "",
    ylimits = (-0.05f0, 1.05f0),
    ax_kwargs = (;),
    legendpos = :rb,
    fig_size = (240, 200),
    fig_theme = CUSTOM_THEME_SMALL,
)
    global solver_props

    if isnothing(solver_fpath_cnames)
        DF.metadata(perf_df, "solver_fpath_cnames", String[])
    end

    with_theme(fig_theme) do
        ## setup figure
        fig = Figure(; size = fig_size)
        ax = Axis(
            fig[1,1];
            title,
            xlabel = L"\tau",
            ylabel = L"\rho",
            limits = (nothing, ylimits),
            yticks = 0f0:0.2f0:1f0,
            ax_kwargs...
        )

        ## filter out `Inf` or `NaN` values
        I = .!( isinf.(perf_df.tau) .|| isnan.(perf_df.tau) )
        tau = perf_df[I, :tau]

        for scname in solver_fpath_cnames 
            sprops = solver_props[scname]
            scatterlines!(
                ax, tau, perf_df[I, scname];
                color = sprops.color,
                label = sprops.name,
                marker = sprops.marker,
                strokecolor = dark_color(sprops.color)
            )
        end
        axislegend(
            ax; 
            position=legendpos,
            font="Latin Modern Mono Light"
        )
        fig
    end
end

# This is how to do a single plot:
let 
    solver_fpath_pair = first(PAIRS_CONSTR)
    purity_df = purity_for_solver_pair(lib_constr, solver_fpath_pair; suffix = "_constr")
    perf_df = MOB.performance_profile(purity_df, solver_fpath_pair)
    plot_performance_profile(perf_df, solver_fpath_pair)
end

# But we have many pairs. 
# And we want to save figures.
# Here a helper: 
function save_fig(fname, fig)
    global PLOTS_PATH
    fig_fpath = joinpath(PLOTS_PATH, fname)
    save(
        ensure_path(fig_fpath), fig;
        px_per_unit = 5, pt_per_unit = 1
    )
    return fig_fpath
end
# And the outer plotting function:
function make_purity_plots(
    results_df, solver_pairs;
    ref_front_cname = "front_fpath", 
    prefix = "", suffix = "",
    title_prefix = "Purity"
)
    global solver_props

    for solver_fpath_pair in solver_pairs
        ## compute purity
        purity_df = purity_for_solver_pair(
            results_df, solver_fpath_pair; 
            ref_front_cname, prefix, suffix
        )

        ## and performance
        perf_df = MOB.performance_profile(purity_df, solver_fpath_pair)
 
        ## generate Figure title
        _sname = sfpath -> haskey(solver_props, sfpath) ? solver_props[sfpath].name : ""
        #src title = "$(title_prefix) – $(_sname(solver_fpath_pair[1])) – $(_sname(solver_fpath_pair[2]))"
        title = title_prefix
        
        ## determine figure save path from purity dataframe metadata
        fig_fpath = joinpath(
            PLOTS_PATH, 
            DF.metadata(purity_df, "results_basename") * ".pdf"
        )
        @info "Making performance plot at `$fig_fpath`."
        perf_fig = plot_performance_profile(
            perf_df, solver_fpath_pair;
            title, fig_size = (160, 150), fig_theme = CUSTOM_THEME_EXTRA_SMALL
        )

        save_fig(fig_fpath, perf_fig)
    end
end

# Finally, make all the purity plots:
make_purity_plots(lib_unconstr, PAIRS_UNCONSTR; 
    suffix = "_unconstr", title_prefix = "Purity (u)")
make_purity_plots(lib_constr, PAIRS_CONSTR; 
    suffix = "_constr", title_prefix = "Purity (c)")
#%%#src
# ## Hypervolume Performance Indicator
# 
# Given some upper limit vector ``r`` (a proxy for the Nadir point), the hypervolume
# for a set of non-dominated points in objective space is the volume enclosed 
# by those points and ``r``.
# The hypervolume is the volume dominated by an approximated Pareto-front.
# The larger the volume, the bigger the hypervolume value.
# The reference front has the best/largest hypervolume value.
# Below we compute the **inverse** of the ratio 
# “hypervolume of approximate solution” ÷ “best hypervolume”.
# The ratio is sometimes called “hyperarea ratio”. 
# By taking the inverse, we have a performance measure that has optimal value 1 and
# takes larger values for worse approximations.

## values for unconstrained test problems
hv_df_uc = MOB.inv_hyperarea_ratio_with_reference(
    lib_unconstr, SOLVERS_UNCONSTR, "front_fpath";
    results_name = "hv_unconstr.jld2",
    viol_tol, 
    overwrite = false,
)
## values for constrained test problems
hv_df_c = MOB.inv_hyperarea_ratio_with_reference(
    lib_constr, SOLVERS_CONSTR, "front_fpath";
    results_name = "hv_constr.jld2",
    viol_tol, 
    overwrite = false,
)
# ”Extract“ the performance profile:
perf_hv_uc = MOB.performance_profile(hv_df_uc, SOLVERS_UNCONSTR)
perf_hv_c = MOB.performance_profile(hv_df_c, SOLVERS_CONSTR)

fig_hv_uc = plot_performance_profile(
    perf_hv_uc, SOLVERS_UNCONSTR;
    title = "Inv. Hyperarea Ratio (u)",
    ax_kwargs = (; limits = ((0.9, 9.99), nothing)) # exclude very large values
)
save_fig("inv_hyperarea_ratio_unconstr.pdf", fig_hv_uc)

fig_hv_c = plot_performance_profile(
    perf_hv_c, SOLVERS_CONSTR;
    title = "Inv. Hyperarea Ratio (c)",
    ax_kwargs = (; limits = ((0.9, 9.99), nothing)) # exclude very large values
)
save_fig("inv_hyperarea_ratio_constr.pdf", fig_hv_c)
#%%#src
# ## Spread Performance
# We can compute different kinds of spread values.
#
# Two possible spread metrics are given in the article by Custódio et al. [^1].
#
# The Ξ (“Xi”) metric measures “spread” by generalizing the 
# two-objective Γ (“Gamma”) metric and giving an idea about 
# holes in the approximated pareto front.
# Ξ is the largest value among distances of consecutive solution value vectors 
# after projecting on the individual objectives.
#
# Θ (“Theta”) generalizes the two-objective Δ (Delta) metric and indicates how much 
# the projected objective distance values vary.
#
# DataFrames for both metrics are created with a single function call 
# to avoid some redundant computations:
## unconstrained problems
Xi_df_uc, Theta_df_uc = MOB.spreads_with_reference(
    lib_unconstr, SOLVERS_UNCONSTR, "front_fpath";
    #lib_unconstr[lib_unconstr.problem_name .== "DPAM1", :], SOLVERS_UNCONSTR, "front_fpath";
    results_name = "spreads_unconstr.jld2",
    viol_tol, 
    overwrite = false,
    normalize_Xi = true
);
## constrained problems
Xi_df_c, Theta_df_c = MOB.spreads_with_reference(
    lib_constr, SOLVERS_CONSTR, "front_fpath";
    results_name = "spreads_constr.jld2",
    viol_tol,
    overwrite = false,
    normalize_Xi = true
)

# Then there is also the "generalized spread” metric.
# It is similar to θ:
Delta_df_uc = MOB.gspread_with_reference(
    lib_unconstr, SOLVERS_UNCONSTR, "front_fpath";
    results_name = "gspread_unconstr.jld2",
    viol_tol, 
    overwrite = false,
)
Delta_df_c = MOB.gspread_with_reference(
    lib_constr, SOLVERS_CONSTR, "front_fpath";
    results_name = "gspread_constr.jld2",
    viol_tol, 
    overwrite = false,
)

# For plotting, we need to turn the value dataframes into 
# performance profile data:
perf_Xi_uc = MOB.performance_profile(Xi_df_uc, SOLVERS_UNCONSTR)
perf_Theta_uc = MOB.performance_profile(Theta_df_uc, SOLVERS_UNCONSTR)
perf_Delta_uc = MOB.performance_profile(Delta_df_uc, SOLVERS_UNCONSTR)

perf_Xi_c = MOB.performance_profile(Xi_df_c, SOLVERS_CONSTR)
perf_Theta_c = MOB.performance_profile(Theta_df_c, SOLVERS_CONSTR)
perf_Delta_c = MOB.performance_profile(Delta_df_c, SOLVERS_CONSTR)

# Finally, make and save the plots:
## unconstrained problems
fig_Theta_uc = plot_performance_profile(perf_Theta_uc, SOLVERS_UNCONSTR; 
    title="Θ Spread (u)", legendpos=(.9, .01))
fig_Xi_uc = plot_performance_profile(perf_Xi_uc, SOLVERS_UNCONSTR; 
    title="Ξ Spread (u)", legendpos=(.01, .4))
fig_Delta_uc = plot_performance_profile(perf_Delta_uc, SOLVERS_UNCONSTR; 
    title="Δ Spread (u)", legendpos=(.9, .01))

save_fig("theta_spread_unconstr.pdf", fig_Theta_uc)
save_fig("xi_spread_unconstr.pdf", fig_Xi_uc)
save_fig("delta_spread_unconstr.pdf", fig_Delta_uc)

## constrained problems
fig_Theta_c = plot_performance_profile(perf_Theta_c, SOLVERS_CONSTR; 
    title="Θ Spread (c)", legendpos=:cb)
fig_Xi_c = plot_performance_profile(perf_Xi_c, SOLVERS_CONSTR; 
    title="Ξ Spread (c)", legendpos=:cb)
fig_Delta_c = plot_performance_profile(perf_Delta_c, SOLVERS_CONSTR; 
    title="Δ Spread (c)", legendpos=:cb)

save_fig("theta_spread_constr.pdf", fig_Theta_c)
save_fig("xi_spread_constr.pdf", fig_Xi_c)
save_fig("delta_spread_constr.pdf", fig_Delta_c)

#%%#src
# ## Constraint Violation

# For problems with nonlinear constraints, we can have a look at the 
# constraint violation values.
# The keywoard argument `quantile_val` is passed to `Statistics.quantile`.
# With `quantile_val = .5` we obtain the median:
mv_c = MOB.quantile_viol(lib_constr, SOLVERS_CONSTR; 
    quantile_val=.5,
)
perf_mv_c = MOB.performance_profile(mv_c, SOLVERS_CONSTR)

fig_mv_c = plot_performance_profile(perf_mv_c, SOLVERS_CONSTR;
    title="Constraint Violation — Median"
)

save_fig("median_constr_viol.pdf", fig_mv_c)

# Now do the same for the third quartile
qv_c = MOB.quantile_viol(lib_constr, SOLVERS_CONSTR; 
    quantile_val=.75,
)
perf_qv_c = MOB.performance_profile(qv_c, SOLVERS_CONSTR)

fig_qv_c = plot_performance_profile(perf_qv_c, SOLVERS_CONSTR;
    title="Constr. Viol. — 3rd Quartile"
)

save_fig("q3_constr_viol.pdf", fig_qv_c)
# ## References
# [^1]: A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
#       “Direct Multisearch for Multiobjective Optimization,” 
#       SIAM J. Optim., vol. 21, no. 3, pp. 1109–1140, Jul. 2011, doi: 10.1137/10079731X.