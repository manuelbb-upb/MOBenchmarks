# ## Helper Callback
# The optimizer allows for a user defined callback.
# Such a callback can specialize `Compromise.check_stopping_criterion`
# to stop the optimization, but also to just gather information.
# In our case, we simply store the current iteration vector at the 
# beginning of each iteration.
# We don't need this operation to be particularly fast, so a vector of vectors
# will do.
# In other experiments, we are also interested in the iteration type, 
# so we store that too:
Base.@kwdef struct InfoGatheringCallback <: C.AbstractStoppingCriterion
    x :: Vector{Vector{Float64}} = Vector{Float64}[]
    it_status :: Vector{C.IterationStatus{Float64}} = C.IterationStatus[]
end

# As indicated, we simply push a copy of the (unscaled) iteration vector
# to the vector of vectors:
function _fill_callback!(cback, vals, iteration_status)
    push!(cback.x, copy(C.cached_ξ(vals)))
    push!(cback.it_status, deepcopy(iteration_status))
    return nothing
end

# In first iteration, store at beginning to get `x0`: 
function C.check_stopping_criterion(
    crit::InfoGatheringCallback,
    stop_point::C.CheckPreIteration,
    mop, mod, scaler, lin_cons, scaled_cons, vals, vals_tmp,
    mod_vals, filter, step_vals, step_cache, crit_cache, trial_caches, 
    iteration_status, iteration_scalars, stop_crits,
    algo_opts
)
    if iteration_scalars.it_index <= 1
        _fill_callback!(crit, vals, iteration_status)
    end
    return nothing
end
# Otherwise, store after iteration:
function C.check_stopping_criterion(
    crit::InfoGatheringCallback,
    stop_point::C.CheckPostIteration,
    mop, mod, scaler, lin_cons, scaled_cons, vals, vals_tmp,
    mod_vals, filter, step_vals, step_cache, crit_cache, trial_caches, 
    iteration_status, iteration_scalars, stop_crits,
    algo_opts
)
    _fill_callback!(crit, vals, iteration_status)
    return nothing
end