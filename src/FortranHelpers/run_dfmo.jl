function run_dfmo!(
    problem_lib_df, 
    dfmo_path
    ;
    kwargs...
) 
    if !ispath(dfmo_path)
        dfmo_path = DF.metadata(problem_lib_df, "dfmo_path", nothing)
    end
    @assert ispath(dfmo_path) "`dfmo_path` not valid."

    df_tmp = run_dfmo(problem_lib_df, dfmo_path; kwargs...)

    colname = DF.metadata(df_tmp, "colname")
    if colname in names(problem_lib_df)
        DF.select!(problem_lib_df, DF.Not(colname))
    end
    DF.leftjoin!(problem_lib_df, df_tmp, on=["problem_name", "constraint_index"])
    DF.metadata!(problem_lib_df, "dfmo_path", dfmo_path; style=:note)

    return problem_lib_df
end

function run_dfmo(
    problem_lib_df, 
    dfmo_path;
    colname = DFMO_RESULTS_COLNAME, 
    log_level=Info,
    indent=0,
    kwargs...
)
    istr = indent_str(indent)

    df_tmp = DF.DataFrame(
        "problem_name" => String[],
        "constraint_index" => Int[],
        colname => String[]
    )

    @progress for row in eachrow(problem_lib_df)
        @unpack problem_name, constraint_index, f90_fpath, so_fpath = row
        @logmsg log_level "$(istr)Querying DFMO results for $(problem_name) / $(constraint_index)."
        results_fpath = run_dfmo(
            f90_fpath, so_fpath, dfmo_path; 
            log_level, indent=indent+1, kwargs...
        )
        isnothing(results_fpath) && continue 
        !isfile(results_fpath) && continue 
        
        push!(df_tmp, (problem_name, constraint_index, results_fpath))
    end
    DF.metadata!(df_tmp, "dfmo_path", dfmo_path; style=:note)
    DF.metadata!(df_tmp, "colname", colname; style=:note)
    return df_tmp
end

function run_dfmo(
    f90_fpath,
    so_fpath,
    dfmo_path
    ;
    results_name = DFMO_RESULTS_FILENAME,
    settings = DFMOSettings(),
    check_settings = true,
    dfmo_exec_filename = DFMO_EXEC_FILENAME,
    overwrite_results = false,
    overwrite_exec = overwrite_results,
    log_level = Info,
    indent = 0
)
    istr = indent_str(indent)

    ## get path to directory with test problem information
    problem_subdir = first(splitdir(f90_fpath))

    ## path for results
    results_fpath = joinpath(problem_subdir, results_name)
    if overwrite_results
        rm(results_fpath; force = true)
    end
        
    if has_valid_field(results_fpath, "settings", settings; check_data = check_settings)
        ## `results_fpath` file exists and stores settings that are 
        ## equal to function argument `settings`
        ## => there is nothing to do and we return early
        return results_fpath
    end

    exec_fpath = joinpath(problem_subdir, dfmo_exec_filename)
    @logmsg log_level "$(istr)Compiling DFMO executable to `$(exec_fpath)`."
    compile_dfmo(
        f90_fpath, exec_fpath, dfmo_path, settings;
        force_recompile = overwrite_exec,
    )
    
    !isfile(exec_fpath) && return nothing
    
    dfmo_dir = joinpath(problem_subdir, "dfmo_results")
    !isdir(dfmo_dir) && mkpath(dfmo_dir)
    @logmsg log_level "$(istr)Running executable."    
    run(Cmd(`$(exec_fpath)`; dir=dfmo_dir), devnull, devnull; wait=true)
    
    @logmsg log_level "$(istr)Remapping `FX`." 
    x, _fx, _viol, num_evals = read_dfmo_results(dfmo_dir)
    results = values_for_matrix(so_fpath, x)            
    JLD2.jldsave(results_fpath; x, num_evals, settings, results...)

    return results_fpath
end

function has_valid_field(
    results_fpath, dataname = "settings", datatarget = settings;
    check_data = true
)
    !isfile(results_fpath) && return false

    !check_data && return true

    old_results = JLD2.load(results_fpath)
    if !haskey(old_results, dataname)
        return false
    end
    
    datatest = old_results[dataname]
    if !isequal(datatarget, datatest)
        return false
    end
    return true
end

function read_dfmo_results(base_path::AbstractString; delete_files=false)
    fobs_path = joinpath(base_path, "pareto_fobs.out")
    vars_path = joinpath(base_path, "pareto_vars.out")
    vars = Vector{Vector{Float64}}()
    fobs = Vector{Vector{Float64}}()
    viol = Vector{Float64}()
    for line in Iterators.drop(eachline(fobs_path), 1)
        vals = [parse(Float64, s) for s in split(line, " ") if !isempty(s)]
        push!(fobs, vals[1:end-1])
        push!(viol, vals[end])
    end
    for line in Iterators.drop(eachline(vars_path), 1)
        vals = [parse(Float64, s) for s in split(line, " ") if !isempty(s)]
        push!(vars, vals[1:end-1])
    end
    F = reduce(hcat, fobs)
    X = reduce(hcat, vars)
   
    fort_path = joinpath(base_path, "fort.1")
    reg =  r"number of function evaluations[^\d]*?(\d+)"
    num_evals = parse(Int, only(match(reg, read(fort_path, String)).captures))

	if delete_files
		rm(fobs_path, force=true)
		rm(vars_path, force=true)
		rm(fort_path, force=true)
	end
    return X, F, viol, num_evals
end


@doc """
    values_for_matrix(evaluator, X; eps=0.1)

Given `evaluator::ProblemEvaluator` and a matrix of input vectors (as columns),
return a `NamedTuple` with fields `x`, `fx`, `cx`, `bx`, `fx_eps` and `viol`.
The fields have the following meaning:
* `x` is a copy of the input matrix
* `fx` is a matrix of objective function values for `evaluator`
* `cx` is a matrix of nonlinear constraint function values for `evaluator`
* `bx` is a matrix of box-constraint adherence values
* `fx_eps` is a matrix of augmented (penalized) objective function values
* `viol` is the vector of constraint violation values

The keyword-argument `eps` specifies the penalty-parameter for the values in 
`fx_eps` and `viol`.

**Only call from within `fortran_wrap`!**
"""
function values_for_matrix(evaluator, X; eps=0.1)
    _n_vars, n_x = size(X)
    @unpack n_vars, n_objfs, n_constrs, lb, ub = evaluator
    @assert _n_vars == n_vars
    
    fx = zeros(n_objfs, n_x)
    cx = zeros(n_constrs, n_x)
    bx = zeros(n_vars, n_x)
    fx_eps = copy(fx)
    viol = zeros(n_x)

    objfs! = inplace_objectives(evaluator)
    constrs! = inplace_constraints(evaluator)
    
    bounds! = function(b, x)
        # x - ub ≤ 0
        # lb - x ≤ 0
        b .= max.( x .- ub, lb .- x )
        nothing
    end
    for (i, x) in enumerate(eachcol(X))
        if any(isnan.(x))
            fill!(@view(fx[:, i]), NaN)
            fill!(@view(cx[:, i]), NaN)
            fill!(@view(bx[:, i]), NaN)
        else
            objfs!(@view(fx[:, i]), x)
            constrs!(@view(cx[:, i]), x)
            bounds!(@view(bx[:, i]), x)
        end
    end
    
    for (i, (y, c, b)) in enumerate(zip(eachcol(fx), eachcol(cx), eachcol(bx)))
        pen = _eps_term(c, eps)
        pen += _eps_term(b, eps)
        @views fx_eps[:, i] .= y .+ pen
        viol[i] = max(0, maximum(c; init=0), maximum(b; init=0))
    end
    return (; x=copy(X), fx, cx, bx, fx_eps, viol)
end

function values_for_matrix(row::DF.DataFrameRow, x)
    return values_for_matrix(row["so_fpath"], x)
end
function values_for_matrix(problem_so::String, x)
    fortran_wrap(problem_so) do evaluator
        values_for_matrix(evaluator, x)
    end
end