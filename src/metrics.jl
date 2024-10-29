REF_FRONT_FILE_VERSION = 0 :: Int
PURITY_RESULT_FILE_VERSION = 0 :: Int
SPREADS_RESULT_FILE_VERSION = 0 :: Int
GSPREAD_RESULT_FILE_VERSION = 0 :: Int
INV_HYPERAREA_RATIO_RESULT_FILE_VERSION = 0 :: Int
QUANTILE_VIOL_FILE_VERSION = 0 :: Int

function performance_profile(
    df_metric_values, solver_fpath_cnames;
)
    solver_fpath_cnames = collect(solver_fpath_cnames)
    num_probs = size(df_metric_values, 1)
    value_matrix = Matrix(df_metric_values[!, solver_fpath_cnames])
    τs = unique(value_matrix[:])
    sort!(τs)
    df_profile = DF.DataFrame("tau" => Float64[], (cname => Float64[] for cname in solver_fpath_cnames)...)
    for τ in τs
        solver_performances = sum.(map(_flag_func(τ), eachcol(value_matrix))) ./ num_probs
        push!(df_profile, [τ; solver_performances])
    end
    _copy_metadata!(df_profile, df_metric_values)
    DF.metadata!(df_profile, "solver_fpath_cnames", solver_fpath_cnames; style=:note)
    return df_profile
end

function _flag_func(τ)
    if isnan(τ)
        return function (svals)
            isnan.(svals)
        end
    end
    return function (svals)
        svals .<= τ
    end
end

function _copy_metadata!(df_trgt, df_src, metafield = "problem_lib_path", metadefault=nothing)
    DF.metadata!(df_trgt, metafield, DF.metadata(df_src, metafield, metadefault); style=:note)
    nothing
end

# # Generalized Spread

# Reference: A. Zhou, Y. Jin, Q. Zhang, B. Sendhoff, and E. Tsang
# Combining model-based and genetics-based offspring generation for
# multi-objective optimization using a convergence criterion,
# 2006 IEEE Congress on Evolutionary Computation, 2006, pp. 3234-3241.

# Note: The reference differs from what is implemented in jmetal and 
# from the description in 
# C. Audet, J. Bigeon, D. Cartier, S. Le Digabel, and L. Salomon, 
# “Performance indicators in multiobjective optimization,” 
# European Journal of Operational Research, vol. 292, no. 2, pp. 397–422, 
# Nov. 2020, doi: 10.1016/j.ejor.2020.11.016.
function gspread_with_reference(
    df, solver_fpath_cnames, ref_front_cname;
    fx_name = "fx",
    viol_name = "viol",
    viol_tol = 1e-3,
    results_path = nothing,
    results_name = "gspread_values",
    overwrite :: Bool = false,
    check_result_file :: Bool = true,
    jmetal_version :: Bool = true
)
    global GSPREAD_RESULT_FILE_VERSION
    res_fpath = _metric_fpath(df, results_path, results_name)
    solver_fpath_cnames = collect(solver_fpath_cnames)
    @assert ref_front_cname in names(df) "Column \"$(ref_front_cname)\" not found."
    if !(is_valid_gspread_file(
        res_fpath, df, solver_fpath_cnames, ref_front_cname; 
        check_result_file, overwrite, jmetal_version
    ))
        hash_df = _empty_hash_df(ref_front_cname, solver_fpath_cnames)   
        num_solvers = length(solver_fpath_cnames)
        Delta_matrix = fill(Inf, num_solvers, size(df, 1))
        @progress for (row_index, row) in enumerate(eachrow(df))
            @unpack problem_name, constraint_index, num_objectives = row
            
            ## initialize extremum variables for performance metric normalization
            Delta_min = Inf
            Delta_max = -Inf

            ## prepare hash variables (to store information for file validity check)
            ref_front_fpath = row[ref_front_cname]
            hash_ref_front = hash_file_contents(ref_front_fpath)
            hash_vec = String[] # stores solver file hashes
            
            ## load reference Pareto front
            F_opt = JLD2.load(ref_front_fpath, fx_name)
            
            ## precompute extremal vectors of reference front
            ## (they are the same for all solvers)
            E = extremal_matrix(F_opt)
       
            ## iterate solver columns
            for (si, scname) in enumerate(solver_fpath_cnames)
                solver_fpath = row[scname]
                
                ## compute and store hash
                push!(hash_vec, hash_file_contents(solver_fpath))

                ## load solver front and filter infeasible points
                F_sol, viol = JLD2.load(solver_fpath, fx_name, viol_name)
                F_sol = F_sol[:, viol .<= viol_tol]

                ## compute spread values and store them
                Delta_matrix[si, row_index] = Delta = generalized_spread(F_opt, F_sol, E; jmetal_version)

                ## update extremum variables
                ## (`Delta` can be `NaN`, which will always evaluate any comparison as `false`)
                if Delta < Delta_min
                    Delta_min = Delta
                end
                if Delta > Delta_max
                    Delta_max = Delta
                end
            end
            if !isinf(Delta_min)
                Delta_matrix[:, row_index] .-= Delta_min
                if Delta_max > Delta_min
                    Delta_matrix[:, row_index] ./= (Delta_max - Delta_min)
                end
            end
            push!(hash_df, [row["problem_name"]; row["constraint_index"]; hash_ref_front; hash_vec])
        end
        df_Delta = _build_res_df(Delta_matrix, solver_fpath_cnames, df)
        if !isa(res_fpath, AbstractString)
            @warn "Cannot save spreads result file. Path is $(res_fpath)."
        else
            JLD2.jldsave(res_fpath;
                file_ver = GSPREAD_RESULT_FILE_VERSION,
                Delta = df_Delta,
                hash_df = hash_df,
                jmetal_version
            )
        end
    else
        df_Delta = JLD2.load(res_fpath, "Delta")
    end

    _copy_metadata!(df_Delta, df)
    DF.metadata!(df_Delta, "fpath", res_fpath; style=:note)
    DF.metadata!(df_Delta, "solver_fpath_cnames", solver_fpath_cnames; style=:note)

    return df_Delta
end


function generalized_spread(
    F_opt :: AbstractMatrix{F}, F_sol :: AbstractMatrix{T}, 
    E = nothing;
    jmetal_version = true
) where {F, T}
    S = Base.promote_type(F, T)
    SNaN = S(NaN)
    SInf = S(Inf)

    _n_objfs, n_opt = size(F_opt)
    if n_opt == 0
        ## if the reference front is empty we can not return something sensible
        return SNaN
    end
    n_objfs, n_sol = size(F_sol)
    if n_sol == 0
        ## if the solver front is empty we can not compute the metric
        return SNaN
    end
    @assert _n_objfs == n_objfs "Mismatch in reference and solver front dimension."

    ## check extremal vectors and recompute them if necessary
    E = _check_extremal_matrix(F_opt, E)

    F_comp = jmetal_version ? F_sol : F_opt
    dists = fill(SInf, size(F_comp, 2))
    for (i, Fi) in enumerate(eachcol(F_comp)) 
        min_dist = dists[i]
        for Fj in eachcol(F_sol)
            dij = sum( (Fi .- Fj).^2 )
            if dij > 0 && dij < min_dist
                min_dist = dij
            end
        end
        dists[i] = isinf(min_dist) ? 0 : sqrt(min_dist)
    end

    dE = zero(S)
    for ei = eachcol(E)
        min_dist = SInf
        for Fj = eachcol(F_sol)
            dij = sum( (ei .- Fj).^2 )
            if dij < min_dist
                min_dist = dij
            end
        end
        dE += sqrt(min_dist)
#        push!(dists, sqrt(min_dist))
    end

    dists_sum = sum(dists)
    dists_avg = dists_sum / n_opt 
    Delta = dE + sum( abs.( dists .- dists_avg ) )
    if Delta > 0
        Delta /= (dE + dists_sum)
    end

    return Delta
end

function _check_extremal_matrix(F, E)
    E_valid = isa(E, AbstractMatrix) && !any(isnan.(E)) && !any(isinf.(E))
    if !E_valid
        return extremal_matrix(F)
    else
        return E
    end
end
function extremal_matrix(F::AbstractMatrix{S}) where S
    n_objfs, n_points = size(F)
    if n_points == 0
        return Matrix{S}(undef, n_objfs, 0)
    end
    E = Matrix{S}(undef, n_objfs, n_objfs)

    for i = 1 : n_objfs
        I = sortperm(F[i, :])
        _i = last(I)
        for j = 1 : n_objfs
            E[j, i] = F[j, _i]
        end
    end
    return E
end
function spreads_with_reference(
    df, solver_fpath_cnames, ref_front_cname;
    fx_name = "fx",
    viol_name = "viol",
    viol_tol = 1e-3,
    results_path = nothing,
    results_name = "spreads_values",
    overwrite :: Bool = false,
    check_result_file :: Bool = true,
    normalize_Xi :: Bool = true,
)
    global SPREADS_RESULT_FILE_VERSION
    res_fpath = _metric_fpath(df, results_path, results_name)
    solver_fpath_cnames = collect(solver_fpath_cnames)
    @assert ref_front_cname in names(df) "Column \"$(ref_front_cname)\" not found."
    if !(is_valid_spreads_file(
        res_fpath, df, solver_fpath_cnames, ref_front_cname; 
        check_result_file, overwrite, normalize_Xi
    ))
        hash_df = _empty_hash_df(ref_front_cname, solver_fpath_cnames)   
        num_solvers = length(solver_fpath_cnames)
        Xi_matrix = fill(Inf, num_solvers, size(df, 1))
        Theta_matrix = similar(Xi_matrix)
        @progress for (row_index, row) in enumerate(eachrow(df))
            @unpack problem_name, constraint_index, num_objectives = row
            
            ## initialize extremum variables for performance metric normalization
            Xi_min = Inf
            Xi_max = -Inf
            Theta_min = Inf
            Theta_max = -Inf

            ## prepare hash variables (to store information for file validity check)
            ref_front_fpath = row[ref_front_cname]
            hash_ref_front = hash_file_contents(ref_front_fpath)
            hash_vec = String[] # stores solver file hashes
            
            ## load reference Pareto front
            F_opt, _F_width, _F_max = JLD2.load(ref_front_fpath, fx_name, "width", "max")
            F_width, F_max = normalize_Xi ? (_F_width, _F_max) : (nothing, nothing)
            ## precompute extremal vectors of reference front
            ## (they are the same for all solvers)
            U, N = extremal_vectors(F_opt)
       
            ## iterate solver columns
            for (si, scname) in enumerate(solver_fpath_cnames)
                solver_fpath = row[scname]
                
                ## compute and store hash
                push!(hash_vec, hash_file_contents(solver_fpath))

                ## load solver front and filter infeasible points
                F_sol, viol = JLD2.load(solver_fpath, fx_name, viol_name)
                F_sol = F_sol[:, viol .<= viol_tol]

                ## compute spread values and store them
                Xi, Theta = Xi_and_Theta_values(F_opt, F_sol, U, N; F_width, F_max)
                Xi_matrix[si, row_index] = Xi
                Theta_matrix[si, row_index] = Theta

                ## update extremum variables
                ## (`Xi` or `Theta` can be `NaN`, which will always evaluate any comparison as `false`)
                if Xi < Xi_min
                    Xi_min = Xi
                end
                if Xi > Xi_max
                    Xi_max = Xi
                end
                if Theta < Theta_min
                    Theta_min = Theta
                end
                if Theta > Theta_max
                    Theta_max = Theta
                end
            end
            if !isinf(Xi_min)
                Xi_matrix[:, row_index] .-= Xi_min
                if Xi_max > Xi_min
                    Xi_matrix[:, row_index] ./= (Xi_max - Xi_min)
                end
            end
            if !isinf(Theta_min)
                Theta_matrix[:, row_index] .-= Theta_min
                if Theta_max > Theta_min
                    Theta_matrix[:, row_index] ./= (Theta_max - Theta_min)
                end
            end
            push!(hash_df, [row["problem_name"]; row["constraint_index"]; hash_ref_front; hash_vec])
        end
        df_Xi = _build_res_df(Xi_matrix, solver_fpath_cnames, df)
        df_Theta = _build_res_df(Theta_matrix, solver_fpath_cnames, df)
        if !isa(res_fpath, AbstractString)
            @warn "Cannot save spreads result file. Path is $(res_fpath)."
        else
            JLD2.jldsave(res_fpath;
                file_ver = SPREADS_RESULT_FILE_VERSION,
                Xi = df_Xi,
                Theta = df_Theta,
                hash_df = hash_df,
                normalize_Xi = normalize_Xi
            )
        end
    else
        df_Xi, df_Theta = JLD2.load(res_fpath, "Xi", "Theta")
    end
   
    _copy_metadata!(df_Xi, df)
    _copy_metadata!(df_Theta, df)
    DF.metadata!(df_Xi, "fpath", res_fpath; style=:note)
    DF.metadata!(df_Theta, "fpath", res_fpath; style=:note)
    DF.metadata!(df_Xi, "solver_fpath_cnames", solver_fpath_cnames; style=:note)
    DF.metadata!(df_Theta, "solver_fpath_cnames", solver_fpath_cnames; style=:note)

    return df_Xi, df_Theta
end

function extremal_vectors(F_matrix::AbstractMatrix{F}) where F
    n_objfs, n_points = size(F_matrix)
    if n_points == 0
        return nothing, nothing
    end
    U = Vector{F}(undef, n_objfs)
    N = similar(U)
    for j = 1 : n_objfs
        u, n = extrema(@view(F_matrix[j, :]))
        U[j] = u
        N[j] = n
    end
    return U, N
end

"""
    Xi_and_Theta_values(F_opt, F_sol, U=nothing, N=nothing)

Given reference Pareto-Front `F_opt` and solver matrix `F_sol`, return 
spread values Ξ (“Xi”) and 𝚯 (“Theta”) as defined in [^1].

The vectors `U`(topia) and `N`(adir) have the extremal points 
of matrix `F_opt` along dimension 2.

[^1:] A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,  
      “Direct Multisearch for Multiobjective Optimization,” 
      SIAM J. Optim., vol. 21, no. 3, pp. 1109–1140, Jul. 2011, doi: 10.1137/10079731X."""
function Xi_and_Theta_values(
    F_opt :: AbstractMatrix{F}, F_sol :: AbstractMatrix{T}, 
    U=nothing, N=nothing;
    F_max = nothing, F_width = nothing
) where {F<:AbstractFloat, T<:AbstractFloat}

    S = Base.promote_type(F, T)
    SNaN = S(NaN)
    SInf = S(Inf)

    _n_objfs, n_opt = size(F_opt)
    if n_opt == 0
        ## if the reference front is empty we can not return something sensible
        return (SNaN, SNaN)
    end
    n_objfs, n_sol = size(F_sol)
    if n_sol == 0
        ## if the solver front is empty we can not compute the metric
        return (SNaN, SNaN)
    end
    @assert _n_objfs == n_objfs "Mismatch in reference and solver front dimension."

    ## check extremal vectors and recompute them if necessary
    U, N = _check_extremal_vectors(F_opt, U, N)
    ## initialize loop variables
    Theta = -SInf
    Xi = -SInf
    dists = zeros(S, n_sol - 1) # distances along objective dimensions

    divisors_Xi = if isnothing(F_max) || isnothing(F_width)
        fill(1.0, n_objfs)
    else
        [iszero(F_width[i]) ? iszero(F_max[i]) ? 1.0 : F_max[i] : F_width[i] for i =1:n_objfs]
    end

    for j = 1 : n_objfs
        ## sort row of solver front in increasing order
        Fj = sort(F_sol[j, :])
        
        ## compute projected distances to extremal points
        Uj = U[j]
        Nj = N[j]

        dist0 = abs(Fj[1] - Uj)
        distN = abs(Nj - Fj[end])

        ## we can already update Xi
        denom_Xi = divisors_Xi[j]
        Xi = max( Xi, dist0 / denom_Xi, distN / denom_Xi )

        ## if there are at least 2 solver points, compute inter-point projected distances
        dists_sum = 0
        for i = 1 : n_sol - 1
            di = Fj[i+1] - Fj[i]
            dists[i] = di
            Xi = max(Xi, di / denom_Xi)                # update Xi with intersite distance
            dists_sum += di
        end
        dists_avg = dists_sum
        if n_sol > 0
            dists_avg /= (n_sol - 1)
        end

        ## Xi is done, now finish Theta
        dists_deviation = sum( abs.(dists .- dists_avg); init = 0)        
       
        dist0N = (dist0 + distN) / 2
        _Theta = abs(dist0 - dist0N) + abs(distN - dist0N) + dists_deviation    # tmp variable
        ## If numerator `_Theta` is positive, then so must be the denominator.
        ## In case `n_sol == 1`, we have a single solver point.
        ## If it matches a single extremal point, then `dist0 + distN == 0` and
        ## we keep `_Theta == 0`,  the best possible value.
        if _Theta > 0
            _Theta /= dist0 + distN + dists_sum
        end

        Theta = max( Theta, _Theta )
    end
    return (Xi, Theta)
end

function _check_extremal_vectors(F_opt, U, N)
    n_objfs = size(F_opt, 1)

    extrema_valid = (
        isa(U, AbstractVector) && length(U) == n_objfs &&
        isa(N, AbstractVector) && length(N) == n_objfs
    )

    if extrema_valid
        for (u, n) in zip(U, N)
            if isinf(u) || isinf(n) || isnan(u) || isnan(n)
                extrema_valid = false
                break
            end
        end
    end

    if !extrema_valid
        U, N = extremal_vectors(F_opt)
    end
    return U, N
end

function purity_with_reference(
    df, solver_fpath_cnames, ref_front_cname;
    fx_name = "fx",
    viol_name = "viol",
    viol_tol = 1e-3,
    results_path = nothing,
    results_name = "purity_values",
    overwrite :: Bool = false,
    check_result_file :: Bool = true
)
    global PURITY_RESULT_FILE_VERSION
    solver_fpath_cnames = collect(solver_fpath_cnames)
    res_fpath = _metric_fpath(df, results_path, results_name)
    @assert ref_front_cname in names(df) "Column \"$(ref_front_cname)\" not found."

    if !(is_valid_purity_file(
        res_fpath, df, solver_fpath_cnames, ref_front_cname; 
        check_result_file, overwrite
    ))
        hash_df = _empty_hash_df(ref_front_cname, solver_fpath_cnames)   
        num_solvers = length(solver_fpath_cnames)
        value_matrix = fill(Inf, num_solvers, size(df, 1))  # initializing with Inf implies that problems for which the reference front is empty, the value is bad (Inf)
        @progress for (row_index, row) in enumerate(eachrow(df))
            ref_front_fpath = row[ref_front_cname]
            hash_ref_front = hash_file_contents(ref_front_fpath)
            hash_vec = String[]
            F_opt = JLD2.load(ref_front_fpath, fx_name)
            n_opt = size(F_opt, 2)
                
            val_min = Inf
            for (cindex, cname) in enumerate(solver_fpath_cnames)
                ## for each solver front `FX`,
                ## determine `n_common`, size of intersection with `F`.
                n_common = 0
                solver_fpath = row[cname]
                push!(hash_vec, hash_file_contents(solver_fpath))
                if n_opt > 0
                    FX, viol = JLD2.load(solver_fpath, fx_name, viol_name)
                    for (i, fx) in enumerate(eachcol(FX))
                        viol[i] > viol_tol && continue
                        for y in eachcol(F_opt)
                            if fx ≈ y
                                n_common += 1
                                break
                            end
                        end
                    end
                    # a large ratio `n_common / n_opt` is good.
                    # purity metric is `n_opt / n_common`, so small values are good
                
                    val_solver = n_opt / n_common
                    if val_solver < val_min
                        val_min = val_solver
                    end

                    value_matrix[cindex, row_index] = val_solver
                end
            end
            if n_opt > 0 && !isinf(val_min)
                value_matrix[:, row_index] ./= val_min
            end
            push!(hash_df, [row["problem_name"]; row["constraint_index"]; hash_ref_front; hash_vec])
        end
        
        df_res = _build_res_df(value_matrix, solver_fpath_cnames, df)
       
        if !isa(res_fpath, AbstractString)
            @warn "Cannot save purity result file. Path is $(res_fpath)."
        else
            JLD2.jldsave(res_fpath;
                file_ver = PURITY_RESULT_FILE_VERSION,
                purity = df_res,
                hash_df = hash_df
            )
        end
    else
        df_res = JLD2.load(res_fpath, "purity")
    end
    
    _copy_metadata!(df_res, df)
    DF.metadata!(df_res, "fpath", res_fpath; style=:note)
    DF.metadata!(df_res, "solver_fpath_cnames", solver_fpath_cnames; style=:note)

    return df_res 
end
# # Median Constraint violation
function quantile_viol(
    df, solver_fpath_cnames;
    viol_name = "viol",
    results_path = nothing,
    results_name = "quantile_viol_values",
    overwrite :: Bool = false,
    check_result_file :: Bool = true,
    quantile_val=.5,
    kwargs...
)
    global QUANTILE_VIOL_FILE_VERSION
    res_fpath = _metric_fpath(df, results_path, results_name)
    solver_fpath_cnames = collect(solver_fpath_cnames)
    if !(is_valid_quantile_viol_file(
        res_fpath, df, solver_fpath_cnames; 
        check_result_file, overwrite,
        viol_name, quantile_val 
    ))
        hash_df = _empty_hash_df(solver_fpath_cnames)   
        num_solvers = length(solver_fpath_cnames)
        mv_matrix = fill(Inf, num_solvers, size(df, 1))
        @progress for (row_index, row) in enumerate(eachrow(df))
            @unpack problem_name, constraint_index, num_objectives = row

            hash_vec = String[]
            mv_min = Inf
            mv_max = -Inf
            for (si, scname) in enumerate(solver_fpath_cnames)
                solver_fpath = row[scname]
                viol = JLD2.load(solver_fpath, viol_name)
                mv = isempty(viol) ? Inf : quantile(viol, quantile_val)
                if mv < mv_min
                    mv_min = mv
                end
                if mv > mv_max
                    mv_max = mv
                end
                push!(hash_vec, hash_file_contents(solver_fpath))
                mv_matrix[si, row_index] = mv
            end
            if !isinf(mv_min)
                mv_matrix[:, row_index] .-= mv_min
                if !isinf(mv_max) && mv_max > mv_min
                    mv_matrix[:, row_index] ./= (mv_max - mv_min)
                end
            end
            push!(hash_df, [problem_name; constraint_index; hash_vec])
        end
        df_mv = _build_res_df(mv_matrix, solver_fpath_cnames, df)
        if !isa(res_fpath, AbstractString)
            @warn "Cannot save median constr. violation result file. Path is $(res_fpath)."
        else
            JLD2.jldsave(
                res_fpath;
                viol_name, quantile_val, hash_df,
                file_ver = QUANTILE_VIOL_FILE_VERSION,
                quantile_viol = df_mv,
            )
        end
    else
        df_mv = JLD2.load(res_fpath, "quantile_viol")
    end

    _copy_metadata!(df_mv, df)
    DF.metadata!(df_mv, "fpath", res_fpath; style=:note)
    DF.metadata!(df_mv, "solver_fpath_cnames", solver_fpath_cnames; style=:note)
    
    return df_mv
end

# # Hypervolume
# HV(F_opt) / HV(F_sol)
# smallest possible: 1, smaller = better
function inv_hyperarea_ratio_with_reference(
    df, solver_fpath_cnames, ref_front_cname;
    fx_name = "fx",
    viol_name = "viol",
    viol_tol = 1e-3,
    results_path = nothing,
    results_name = "inv_hyperarea_ratio_values",
    overwrite :: Bool = false,
    check_result_file :: Bool = true,
)
    global INV_HYPERAREA_RATIO_RESULT_FILE_VERSION
    res_fpath = _metric_fpath(df, results_path, results_name)
    solver_fpath_cnames = collect(solver_fpath_cnames)
    @assert ref_front_cname in names(df) "Column \"$(ref_front_cname)\" not found."
    if !(is_valid_inv_hyperarea_ratio_file(
        res_fpath, df, solver_fpath_cnames, ref_front_cname; 
        check_result_file, overwrite, viol_tol, viol_name
    ))
        hash_df = _empty_hash_df(ref_front_cname, solver_fpath_cnames)   
        num_solvers = length(solver_fpath_cnames)
        hv_matrix = fill(Inf, num_solvers, size(df, 1))
        @progress for (row_index, row) in enumerate(eachrow(df))
                      
            ## prepare hash variables (to store information for file validity check)
            ref_front_fpath = row[ref_front_cname]
            hash_ref_front = hash_file_contents(ref_front_fpath)
            hash_vec = String[] # stores solver file hashes
            
            ## load reference Pareto front
            F_opt, ref_vec = JLD2.load(ref_front_fpath, fx_name, "max")
            
            hv_opt = hypervolume(Matrix(F_opt'), ref_vec)

            ## initialize extremum variables for performance metric normalization
            hv_min = hv_opt > 0 ? Inf : 1

            ## iterate solver columns
            for (si, scname) in enumerate(solver_fpath_cnames)
                solver_fpath = row[scname]
                
                ## compute and store hash
                push!(hash_vec, hash_file_contents(solver_fpath))

                hv = 1
                if hv_opt > 0
                    ## load solver front and filter infeasible points
                    F_sol, viol = JLD2.load(solver_fpath, fx_name, viol_name)
                    F_sol = F_sol[:, viol .<= viol_tol]
                    I = pareto_index(F_sol)
                    F_sol = Matrix(F_sol[:, I]')

                    ## compute spread values and store them
                    hv_sol = hypervolume(F_sol, ref_vec)
                    hv = hv_opt / hv_sol
                    ## update extremum variables
                    if hv < hv_min
                        hv_min = hv
                    end
                end

                hv_matrix[si, row_index] = hv
            end
            hv_matrix[:, row_index] ./= hv_min
                
            push!(hash_df, 
                [row["problem_name"]; row["constraint_index"]; hash_ref_front; hash_vec])
        end
        df_hv = _build_res_df(hv_matrix, solver_fpath_cnames, df)
        if !isa(res_fpath, AbstractString)
            @warn "Cannot save HV result file. Path is $(res_fpath)."
        else
            JLD2.jldsave(
                res_fpath;
                viol_name,
                viol_tol,
                file_ver = INV_HYPERAREA_RATIO_RESULT_FILE_VERSION,
                inv_hyperarea_ratio = df_hv,
                hash_df = hash_df,
            )
        end
    else
        df_hv = JLD2.load(res_fpath, "inv_hyperarea_ratio")
    end

    _copy_metadata!(df_hv, df)
    DF.metadata!(df_hv, "fpath", res_fpath; style=:note)
    DF.metadata!(df_hv, "solver_fpath_cnames", solver_fpath_cnames; style=:note)

    return df_hv
end

function _build_res_df(value_matrix, solver_fpath_cnames, df)
    df_res = copy(df)
    
    for (j, s_cname) in enumerate(solver_fpath_cnames)
        sname = string(s_cname)
        DF.rename!(df_res, sname => "_" * sname)
        df_res[:, sname] = vec(value_matrix[j, :])
    end

    DF.select!(
        df_res,
        "problem_name",
        "constraint_index",
        solver_fpath_cnames,
        DF.Not([solver_fpath_cnames; "problem_name"; "constraint_index"])
    )
   
    return df_res
end

function _empty_hash_df(ref_front_cname, solver_fpath_cnames)
    solver_fpath_cnames = [ref_front_cname; solver_fpath_cnames]
    return _empty_hash_df(solver_fpath_cnames)
end
function _empty_hash_df(solver_fpath_cnames)
    return DF.DataFrame(
        "problem_name" => String[],
        "constraint_index" => Int[],
        (cname => String[] for cname in solver_fpath_cnames)...
    )
end

function _metric_fpath(df, results_path, results_name)
    if !isa(results_path, AbstractString)
        results_path = DF.metadata(df, "problem_lib_path", nothing)
    end
    if !isa(results_name, AbstractString)
        results_name = "purity_values.jld2"
    end
    res_fpath = if isa(results_path, AbstractString)
        if !isdir(results_path)
            mkpath(results_path)
        end
        abspath(joinpath(results_path, results_name))
    else
        nothing
    end
    return res_fpath
end

function is_valid_quantile_viol_file(
    pur_res_fpath, df, solver_fpath_cnames; 
    check_result_file, overwrite,
    viol_name, quantile_val,
)
    global QUANTILE_VIOL_FILE_VERSION
    return is_valid_result_file_with_reference(
        pur_res_fpath, df, solver_fpath_cnames, nothing;
        check_result_file, overwrite,
        FILE_VER = QUANTILE_VIOL_FILE_VERSION,
        viol_name, quantile_val
    )
end
function is_valid_inv_hyperarea_ratio_file(
    pur_res_fpath, df, solver_fpath_cnames, ref_front_cname; 
    check_result_file, overwrite, viol_name, viol_tol
)
    global INV_HYPERAREA_RATIO_RESULT_FILE_VERSION
    is_valid = is_valid_result_file_with_reference(
        pur_res_fpath, df, solver_fpath_cnames, ref_front_cname;
        check_result_file, overwrite, 
        FILE_VER = INV_HYPERAREA_RATIO_RESULT_FILE_VERSION,
        viol_name, viol_tol
    )
    return is_valid
end
function is_valid_gspread_file(
    res_fpath, df, solver_fpath_cnames, ref_front_cname; 
    check_result_file, overwrite,
    jmetal_version
)
    global GSPREAD_RESULT_FILE_VERSION 
    is_valid = is_valid_result_file_with_reference(
        res_fpath, df, solver_fpath_cnames, ref_front_cname;
        check_result_file, overwrite, 
        FILE_VER = GSPREAD_RESULT_FILE_VERSION
    )
    !is_valid && return is_valid
    return jmetal_version == JLD2.load(res_fpath, "jmetal_version")
end

function is_valid_spreads_file(
    res_fpath, df, solver_fpath_cnames, ref_front_cname; 
    check_result_file, overwrite,
    normalize_Xi
)
    global SPREADS_RESULT_FILE_VERSION 
    is_valid = is_valid_result_file_with_reference(
        res_fpath, df, solver_fpath_cnames, ref_front_cname;
        check_result_file, overwrite, 
        FILE_VER = SPREADS_RESULT_FILE_VERSION
    )
    !is_valid && return is_valid
    return normalize_Xi == JLD2.load(res_fpath, "normalize_Xi")
end
function is_valid_purity_file(
    pur_res_fpath, df, solver_fpath_cnames, ref_front_cname; 
    check_result_file, overwrite
)
    global PURITY_RESULT_FILE_VERSION
    return is_valid_result_file_with_reference(
        pur_res_fpath, df, solver_fpath_cnames, ref_front_cname;
        check_result_file, overwrite, 
        FILE_VER = PURITY_RESULT_FILE_VERSION
    )
end

function is_valid_result_file_with_reference(
    res_fpath, df, solver_fpath_cnames, ref_front_cname = nothing; 
    check_result_file, overwrite, FILE_VER,
    kwargs...
)
    !isa(res_fpath, AbstractString) && return false
    if overwrite
        rm(res_fpath; force=true)
    end
    if !isfile(res_fpath)
        return false
    end
    
    if check_result_file
        hashed_cnames = isnothing(ref_front_cname) ? solver_fpath_cnames : vcat(solver_fpath_cnames, ref_front_cname)
        file_ver, hash_df = JLD2.load(res_fpath, "file_ver", "hash_df")
        if file_ver != FILE_VER
            return false
        end
        if size(hash_df, 1) != size(df, 1)
            return false
        end

        kw_trgts = JLD2.load(res_fpath, string.(keys(kwargs))... )
        for (i, kw_val) in enumerate(values(kwargs))
            if kw_val != kw_trgts[i]
                return false
            end
        end
        I = sortperm(df, ["problem_name", "constraint_index"])
        I = invperm(I)
        DF.sort!(hash_df, ["problem_name", "constraint_index"])
        hash_df = hash_df[I, :]
        for cname in hashed_cnames
            if !(cname in names(hash_df))
                return false
            end
            ref_hash = hash_df[!, cname]
            is_hash = map(
                fpath -> hash_file_contents(fpath),
                df[!, cname]
            )
            if ref_hash != is_hash
                return false
            end
        end
    end
    return true
end

function combine_fronts!(
    df::DF.DataFrame,
    solver_fpath_cnames,
    ;
    cname ::AbstractString = "front_fpath",
    kwargs...
)
   df_tmp = combine_fronts(df, solver_fpath_cnames; cname, kwargs...)

   if cname in names(df)
    DF.select!(df, DF.Not(cname))
   end
   DF.leftjoin!(df, df_tmp; on = ["problem_name", "constraint_index"])
   return df
end

function combine_fronts(
    df::DF.DataFrame,
    solver_fpath_cnames;
    cname ::AbstractString = "front_fpath",
    kwargs...
)
    @assert "f90_fpath" in names(df)
    df_tmp = DF.DataFrame(
        "problem_name" => String[],
        "constraint_index" => Int[],
        cname => String[],
    )
    @progress for row = eachrow(df)    
        front_fpath = combine_fronts(row, solver_fpath_cnames; kwargs...)
        @unpack problem_name, constraint_index = row
        push!(df_tmp, (problem_name, constraint_index, front_fpath))
    end
    return df_tmp
end

function combine_fronts(
    row::DF.DataFrameRow,
    solver_fpath_cnames;
    results_name = "combined_front.jld2",
    fx_name = "fx",
    viol_name = "viol",
    viol_tol = 1e-3,
    overwrite = false,
    check_front = true
)
    global REF_FRONT_FILE_VERSION
    @unpack problem_name, constraint_index, num_objectives, f90_fpath = row
        
    subdir = first(splitdir(f90_fpath))
    front_fpath = joinpath(subdir, results_name)

    if !is_valid_front_file(
        front_fpath, row, solver_fpath_cnames;
        check_front, overwrite
    )
        F_all = Matrix{Float64}(undef, num_objectives, 0)
        row_min = fill(Inf, num_objectives)
        row_max = fill(-Inf, num_objectives)
        hash_dict = Dict{String, String}()
        for cname in solver_fpath_cnames
            fpath = row[cname]
            col_hash = hash_file_contents(fpath)
            hash_dict[cname] = col_hash
            _F_col, v_col = JLD2.load(fpath, fx_name, viol_name)

            _rmin, _rmax = let row_extrema = extrema(_F_col; dims=2);
                first.(row_extrema), last.(row_extrema)
            end
            row_min .= min.(row_min, _rmin)
            row_max .= max.(row_max, _rmax)

            F_col = _F_col[:, v_col .<= viol_tol ]
            F_all = hcat(F_all, F_col)
        end
        F = pareto_union(F_all)
        JLD2.jldopen(front_fpath, "w") do jldf
            jldf[fx_name] = F
            jldf["min"] = row_min
            jldf["max"] = row_max
            jldf["width"] = row_max .- row_min
            jldf["file_ver"] = REF_FRONT_FILE_VERSION
            jldf["hash_dict"] = hash_dict
        end
    end
    return front_fpath
end

function is_valid_front_file(
    front_fpath, row, solver_fpath_cnames; 
    check_front, overwrite
)
    global REF_FRONT_FILE_VERSION
    if overwrite 
        rm(front_fpath; force=true)
    end
    if !isfile(front_fpath)
        return false
    end
    
    if check_front
        file_ver, hash_dict = JLD2.load(front_fpath, "file_ver", "hash_dict")
        if file_ver != REF_FRONT_FILE_VERSION
            return false
        end
        for cname in solver_fpath_cnames
            if !haskey(hash_dict, cname)
                return false
            end
            ref_hash = hash_dict[cname]
            is_hash = hash_file_contents(row[cname])
            if ref_hash != is_hash
                return false
            end
        end
    end
    return true
end

function pareto_union(Fs::Vararg{<:AbstractMatrix})
    F = reduce(hcat, Fs)
    isnondom = pareto_index(F)
    return F[:, isnondom]
end

pareto_index(F::AbstractMatrix)=pareto_index(F, nothing)

function pareto_index(
    F::AbstractMatrix, isdom::Nothing
)
    isdom = zeros(Bool, size(F, 2))
    return pareto_index(F, isdom)
end

function pareto_index(
    F::AbstractMatrix, isdom::AbstractVector{Bool}
)
    N = size(F, 2)
    for j = 1:N
        isdom[j] && continue
        fj = @view(F[:, j])
        for i = 1:N
            i == j && continue
            isdom[i] && continue
            fi = @view(F[:, i])
            if fi ⪯ fj
                ## by using `⪯` instead of `≺` we filter duplicates as well
                isdom[j] = true
                break
            end
            if fj ≺ fi
                isdom[i] = true
            end
        end
    end
    return .!(isdom)
end

⪯(a, b) = vecleq(a, b)
vecleq(a, b) = all( a .<= b )

≺(a, b) = dominates(a, b)
function dominates(a, b)
    if a ⪯ b
        return any( a .< b )
    end
    return false
end
function isdominatedby(a, b)
    return dominates(b, a)
end

