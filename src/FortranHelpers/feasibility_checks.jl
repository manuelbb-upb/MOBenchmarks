# # Third Stage -- Feasibility Check(s)

"""
    FeasibilityTestSettings(; <keyword arguments>)

A configuration object for the NLopt feasibility test.

# Keyword arguments

* `stopval::Float64=1e-3` -- Satisfactory constraint violation value.  
  Optimizer stops if this value is reached.
* `tol::Float64=1e-6` -- Tolerance for NLopt inequality constraints.
* `maxtime_factor::Int=2` -- Optimizer runs for `maxtime_factor * N` seconds, where  
  `N` is the number of variables.
* `algorithm::Symbol=:LN_COBYLA` -- NLopt algorithm (has to be derivative-free)."""
@with_kw struct FeasibilityTestSettings
    stopval :: Float64 = 1e-3
    tol :: Float64 = min(1e-3, stopval^2)
    maxtime_factor :: Int = 2
    algorithm :: Symbol = :LN_COBYLA

    @assert string(algorithm)[2] == 'N' "NLopt algorithm should be derivative-free."
end
@batteries FeasibilityTestSettings selfconstructor=false


function check_feasibility!(
    problem_lib_df;
    kwargs...
)
    df_tmp = check_feasibility(problem_lib_df; kwargs...)

    cnames = setdiff(names(df_tmp), ["problem_name", "constraint_index"])
    DF.select!(problem_lib_df, setdiff(names(problem_lib_df), cnames))
    
    DF.leftjoin!(problem_lib_df, df_tmp; on = ["problem_name", "constraint_index"])
    DF.dropmissing!(problem_lib_df)
    return problem_lib_df
end 

function check_feasibility(
    problem_lib_df;
    problem_lib_path=nothing,
    feasibility_test_settings = FeasibilityTestSettings(),
    feasibility_table_subdir = FEASIBILITY_LIB_SUBDIR, 
    feasibility_info_filename = FEASIBILITY_INFO_FILENAME,
    overwrite_feasibility_subdir = false,
    overwrite_feasibility_subdir_files = overwrite_feasibility_subdir,
    overwrite_feasibility_info_files = false,
    cnames = ["feas_fpath",],
    only_feasible = false,
    log_level=Info,
    indent=0,
    kwargs...
)   
    ## initialize 2 dataframes:
    ## `table` has columns "constraint_index", "lb", "ub" and feasibility information columns
    ## `df_tmp` is filled with rows corresponding to `problem_lib_df`
    ## We use `table` to avoid unnecessary checks for equivalent feasibility problems
    table, df_tmp = _init_feasibility_check_dataframes(
        problem_lib_df, problem_lib_path, feasibility_table_subdir, overwrite_feasibility_subdir
    )

    problem_lib_path = DF.metadata!(df_tmp, "problem_lib_path", nothing)
    isnothing(problem_lib_path) && return problem_lib_path
 
    table_lib_path = DF.metadata(table, "feasibility_table_lib_path")
 
    @progress for row in eachrow(problem_lib_df)
        _check_row_feasibility!(
            df_tmp, table,
            row, feasibility_test_settings, table_lib_path, 
            overwrite_feasibility_subdir_files, feasibility_info_filename, 
            overwrite_feasibility_info_files; 
            log_level, indent
        )
       
    end

    if only_feasible
        @unpack stopval = feasibility_test_settings
        df_tmp = df_tmp[ df_tmp.theta .<= stopval, : ]
    end
    DF.select!(df_tmp, "problem_name", "constraint_index", cnames)
    return df_tmp
end

function _init_feasibility_check_dataframes(
    problem_lib_df, 
    problem_lib_path,
    feasibility_table_subdir,
    overwrite_feasibility_subdir
)
    if isnothing(problem_lib_path)
        problem_lib_path = DF.metadata(problem_lib_df, "problem_lib_path", nothing)
    end
    
    table_lib_path = joinpath(problem_lib_path, feasibility_table_subdir) 
    if overwrite_feasibility_subdir
        rm(table_lib_path; force=true, recursive=true)
    end
    if !isdir(table_lib_path)
        mkpath(table_lib_path)
    end
    
    table = _empty_df_feasibility_table()
    df_tmp = _empty_df_feasibility_results()
    DF.metadata!(table, "feasibility_table_lib_path", table_lib_path)
    DF.metadata!(df_tmp, "problem_lib_path", problem_lib_path)
    return table, df_tmp
end

function _empty_df_feasibility_table()
    return DataFrame(
        "constraint_index" => Int[],
        "box_hash" => UInt64[],
        "xfeas" => Vector{Float64}[],
        "theta" => Float64[],
        "nlopt_ret" => Symbol[],
    )
end

function _empty_df_feasibility_results()
    return DataFrame(
        "problem_name" => String[],
        "constraint_index" => Int[],
        "feas_fpath" => String[],
        "x0_feasible" => Bool[],
        "xfeas" => Vector{Float64}[],
        "theta" => Float64[],
        "nlopt_ret" => Symbol[],
    )
end

function _check_row_feasibility!(
    df_tmp, table, 
    row,
    feasibility_test_settings, 
    table_lib_path,
    overwrite_feasibility_subdir_files,
    feasibility_info_filename, 
    overwrite_feasibility_info_files;
    log_level, indent
)
    ## unpack information from row
    @unpack problem_name, constraint_index, num_vars, num_constraints, f90_fpath = row
    x0, lb, ub = _row_dat(row, "x0", "lb", "ub")

    ## log information to user
    istr = indent_str(indent)
    @logmsg log_level "$(istr)Checking feasibility of $(problem_name) / $(constraint_index)."

    ## ensure constraint index subdir exists, e.g., `PROBLEMLIB/FEASIBILITY/0`
    subdir = _feasibility_constraint_subdir(table_lib_path, constraint_index)

    ## ensure subfolder for current box constraint set `(lb, ub)` exists,
    ## and within that subfolder there is a file for current starting point `x0`:
    box_hash, xbox_fpath = _feasiblity_xbox(
        x0, lb, ub, subdir, overwrite_feasibility_subdir_files
    )
   
    ## ensure that in the problem library, in the same folder as `f90_fpath`,
    ## we can put a feasibility info file 
    problem_feasibility_fpath = _feasibility_feas_fpath(
        f90_fpath, feasibility_info_filename, overwrite_feasibility_info_files
    )

    ## the above function call delete information files if requested;
    ## if any of those files still exists, we use the data stored therein
    src_fpath = if isfile(xbox_fpath)
        xbox_fpath
    elseif isfile(problem_feasibility_fpath)
        problem_feasibility_fpath
    else
        nothing
    end

    (x0_feasible, xfeas, theta, nlopt_ret) = if !isnothing(src_fpath)
        ## use data from file
        JLD2.load(
            src_fpath, "x0_feasible", "xfeas", "theta", "nlopt_ret")
    else
        ## we have to check feasibility because there is no data stored
        _feasibility_perform_check!(
            table,
            x0, lb, ub, box_hash, constraint_index, num_constraints, feasibility_test_settings;
            log_level, indent
        )
    end
    push!(
        df_tmp, 
        (; 
            problem_name, constraint_index, feas_fpath = problem_feasibility_fpath,
            x0_feasible, xfeas, theta, nlopt_ret
        )
    )

    if !isfile(xbox_fpath) 
        JLD2.jldsave(xbox_fpath; x0_feasible, xfeas, theta, nlopt_ret)
    end
    ## `xbox_fpath` has priority over `feas_fpath`, and we make sure to keep them in sync
    cp(xbox_fpath, problem_feasibility_fpath; force = true) 
    return nothing
end

function _feasibility_constraint_subdir(table_lib_path, constraint_index)
    subdir = joinpath(table_lib_path, string(constraint_index))
    if !isdir(subdir)
        mkpath(subdir)
    end
    return subdir
end

function _feasibility_feas_fpath(
    f90_fpath, feasibility_info_filename, overwrite_feasibility_info_files
)
    problem_base_dir = first(splitdir(f90_fpath))
    feas_fpath = joinpath(problem_base_dir, feasibility_info_filename)
    if overwrite_feasibility_info_files
        rm(feas_fpath; force=true)
    end
    return feas_fpath
end

function _feasiblity_xbox(
    x0, lb, ub, subdir,
    overwrite_feasibility_subdir_files
)
    box_hash = hash(lb, hash(ub)) 
    box_dir = joinpath(subdir, string(box_hash))
    if !isdir(box_dir)
        mkpath(box_dir)
    end

    xbox_hash = hash(x0, box_hash)
    xbox_fpath = joinpath(box_dir, string(xbox_hash) * ".jld2")

    if overwrite_feasibility_subdir_files
        rm(xbox_fpath; force = true)
    end
    
    return box_hash, xbox_fpath
end

function _feasibility_perform_check!(
    table,  # possibly modified
    x0, lb, ub, box_hash, constraint_index, num_constraints, feasibility_test_settings;
    log_level, indent
)
    @unpack stopval = feasibility_test_settings
    constrs! = julia_constraint_func(Val(constraint_index))
    x0_feasible = _feasibility_pretest(x0, lb, ub, num_constraints, constrs!; stopval)
    if x0_feasible
        xfeas = x0
        theta = 0
        nlopt_ret = :NotNeeded
    else
        i = findfirst( 
            (table.constraint_index .== constraint_index) .& 
            (table.box_hash .== box_hash)
        )
        if isnothing(i)
            @unpack tol, maxtime_factor, algorithm = feasibility_test_settings
            xfeas, theta, nlopt_ret = smooth_feasibility_test(
                x0, lb, ub, num_constraints, constrs!; 
                stopval, tol, maxtime_factor, algorithm,
                log_level, indent = indent + 1
            )
            push!(table, (; box_hash, constraint_index, xfeas, theta, nlopt_ret))
        else
            @unpack xfeas, theta, nlopt_ret = table[i, :]
        end
    end
    return (x0_feasible, xfeas, theta, nlopt_ret)
end

function _feasibility_pretest(
    x0, lb, ub, n_constrs, @nospecialize(constrs!);
    stopval = 1e-3
)
    if n_constrs <= 0
        return all(lb .<= x0 .<= ub)
    end

    if all(lb .<= x0 .<= ub)
        stopval = max(stopval, 0)
        
        gx = zeros(n_constrs)
        constrs!(gx, x0)
        theta = maximum(gx)
        if theta .<= stopval
            return true
        end
    end
    return false
end
        
function smooth_feasibility_test(
    x0, lb, ub, n_constrs, @nospecialize(constrs!); 
    stopval :: Float64 = 1e-3,
    tol :: Float64 = stopval^2, 
    maxtime_factor :: Int = 2,
    algorithm :: Symbol = :LN_COBYLA,
    log_level = Info,
    indent = 0
) 
    istr = indent_str(indent)
    stopval = max(stopval, 0)
        
    gx = zeros(n_constrs)
    constrs!(gx, x0)
    theta = maximum(gx)
    
    n_vars = length(x0)
    _lb = [min(-2*abs(theta), -1); lb]
    _ub = [max(2*abs(theta), 1); ub]
        
    opt = NLopt.Opt(algorithm, n_vars + 1)
    opt.lower_bounds = _lb
    opt.upper_bounds = _ub 
    nlopt_objf = function(tx, grad_vec)
        return first(tx)
    end
    opt.min_objective = nlopt_objf

    nlopt_constraint_func! = function(constr_vec, tx, grads)
        t = first(tx)
        x = tx[2:end]
        constrs!(gx, x)
        constr_vec .= gx .- t
        return constr_vec
    end
    NLopt.inequality_constraint!(opt, nlopt_constraint_func!, fill(tol, n_constrs))

    maxtime = maxtime_factor * n_vars
    opt.maxtime = maxtime
    opt.stopval = stopval
        
    tx0 = [theta; x0]

    @logmsg log_level "$(istr)Starting optimization for at most $(maxtime) secs."
    (optf, opt_tx, ret) = NLopt.optimize(opt, tx0)
    
    xopt = opt_tx[2:end]
    xopt = min.(max.(lb, xopt), ub)
    constrs!(gx, xopt)
    theta = max(maximum(gx), 0)
    
    @logmsg log_level "$(istr)Finished with $(ret) and theta $theta."
    return xopt, theta, ret 
end