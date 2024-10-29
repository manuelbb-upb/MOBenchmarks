# Load packages and definitions
include("preamble.jl")

# Initialize and load CompromiseHelpers extension:
CH = MOB.init_CompromiseHelpers()
import .CH: SettingsContainer, run_compromise!

#%%#src
# Load problem library information from CSV
lib_df = MOB.readlib(PROBLEM_LIB_CSV_FPATH);

# # Experiments

# Make the number of initial points and maximum model points dependent on number of 
# variables by adding columns to dataframe:
max_func_calls = 2000
max_num_halton = ceil(Int, max_func_calls / 4)
lib_df.num_halton = min.( max_num_halton, lib_df.num_vars * 100 )
lib_df.max_points = lib_df.num_vars .* 5
#%%#src
# ## Optimizer Settings
# Specify settings for testing...
# Inner options if there is an outer algorithm:
common_opts = C.AlgorithmOptions(;
    log_level = Debug,
    stop_delta_min = 1e-10,
    #trial_mode = Val(:max_diff),
    trial_mode = Val(:max_rho),
    trial_update = Val(:stepsize),
    nu_accept = 1e-4,
    nu_success = 0.9,
    step_config = C.SteepestDescentConfig(;
        #backtracking_mode = Val(:max),
        backtracking_mode = Val(:any),
        rhs_factor = 0.01,
        normal_step_norm = 1,
        descent_step_norm = 1,
        normalize_gradients = true,
    )
)
# (Outer) Settings for sequential optimization runs:
algo_opts_seq = C.SequentialOuterAlgorithmOptions(;
    log_level = Debug,
    sort_by_delta = true,                   # try points with large radius first
    delta_factor = 0.05,                    # skip solutions with small radius
    initial_nondominance_testing = false,   # only consider nondominated points
    nondominance_testing_offset = 10,       # filter points in some iterations or not
    final_nondominance_testing = true,      # return only nondominated points
    inner_opts = common_opts
)
## wrap in a `SettingsContainer` to enable dispatch without type piracy:
compromise_sequential_settings = SettingsContainer(;
    num_halton = "num_halton",
    max_points = "max_points",
    max_func_calls,
    kernel = C.CubicKernel(),
    eps = nothing,
    algo_opts = algo_opts_seq
)
# Use the same outer settings for penalized sequential optimization with nonlinear constraints:
compromise_penalized_settings = SettingsContainer(;
    num_halton = "num_halton",
    max_points = "max_points",
    max_func_calls,
    kernel = C.CubicKernel(),
    eps = 0.1,
    algo_opts = algo_opts_seq
)

# Now for the set-based approach:
compromise_set_settings = SettingsContainer(;
    num_halton = "num_halton",
    max_points = "max_points",
    max_func_calls,
    kernel = C.CubicKernel(),
    eps = nothing,
    algo_opts = common_opts
)
#%%#src
# ## Meta Settings
overwrite_results = false
disable_logging = false
thread_rows = true
log_level=Info
inner_log_level=Debug

#%%#src
# ## Start Optimization
# Run set-based optimization for all chosen test problems:
@info "Running set-based optimizer ..."
run_compromise!(
    lib_df,
    ;
    colname=COL_NAME_COMPROMISE_SET,
    results_name="cubic_cset.jld2",
    logging_prefix="SET ", 
    settings=compromise_set_settings,
    overwrite_results,
    disable_logging,
    thread_rows,
    log_level, 
    inner_log_level
)
#%%#src
# Run sequential optimization for all chosen test problems:
@info "Running sequential optimizer ..."
run_compromise!(
    lib_df;
    colname = COL_NAME_COMPROMISE_SEQ,
    results_name = "cubic_cseq.jld2", 
    logging_prefix="SEQ ", 
    settings=compromise_sequential_settings,
    overwrite_results,
    disable_logging,
    thread_rows,
    log_level, 
    inner_log_level
)
#%%#src
# Run set-based optimization for all constrained test problems:
lib_unconstr = lib_df[lib_df.constraint_index .== 0, :]
lib_constr = lib_df[lib_df.constraint_index .> 0, :]
@info "Running augmented optimizer for constrained problems..."
run_compromise!(
    lib_constr,
    ;
    colname=COL_NAME_COMPROMISE_PEN, 
    results_name="cubic_cpen.jld2", 
    logging_prefix="PEN ",
    settings=compromise_penalized_settings,
    overwrite_results,
    disable_logging,
    thread_rows,
    log_level, 
    inner_log_level,
)
# ## Store Results
# The dataframes `lib_constr` and `lib_unconstr` have been augmented with columns pointing
# to result files.
# Save this information:
MOB.store!(lib_constr; outname=LIB_CONSTR_CSV_NAME)
MOB.store!(lib_unconstr; outname=LIB_UNCONSTR_CSV_NAME)