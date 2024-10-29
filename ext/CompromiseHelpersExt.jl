module CompromiseHelpersExt

import MOBenchmarks as MOB
import MOBenchmarks.FortranHelpers as BFH
import MOBenchmarks.FortranHelpers: fortran_wrap, inplace_objectives, 
    inplace_constraints, l1_penalized_objectives, lu_penalized_objectives,
    has_valid_field, values_for_matrix

import Compromise
import Compromise: MutableMOP, AlgorithmOptions, SequentialOuterAlgorithmOptions,
    add_objectives!, add_nl_ineq_constraints!, optimize_with_algo, optimize_many,
    CubicKernel

import LoggingExtras: TransformerLogger
import Accessors: @reset
import DataFrames as DF
import UnPack: @unpack
import StructHelpers: @batteries
import HaltonSequences: HaltonPoint

import JLD2

import Logging
import Logging: @logmsg, Info
import Dates

import ProgressLogging: @progress, @withprogress, @logprogress

Base.@kwdef struct SettingsContainer{
    num_halton_Type <: Union{Integer, Symbol, AbstractString},
    max_func_calls_Type <: Union{Integer, Symbol, AbstractString, Nothing},
    max_points_Type <: Union{Integer, Symbol, AbstractString, Nothing},
    algo_opts_Type,
    kernel_Type,
    eps_Type <: Union{Real, Nothing}
}
    num_halton :: num_halton_Type = 0
    max_func_calls :: max_func_calls_Type = 2000
    max_points :: max_points_Type = nothing
    algo_opts :: algo_opts_Type = SequentialOuterAlgorithmOptions()
    kernel :: kernel_Type = CubicKernel()
    eps :: eps_Type = nothing
end
@batteries SettingsContainer selfconstructor=false

_compromise_name(::Any)="comp_results"
_compromise_name(settings::SettingsContainer)=_compromise_name(settings.algo_opts, settings.eps)
_compromise_name(::Any, ::Any)="comp_results"
_compromise_name(algo_opts::AlgorithmOptions, eps)="comp_set"
_compromise_name(algo_opts::SequentialOuterAlgorithmOptions, eps::Nothing)="comp_seq"
_compromise_name(algo_opts::SequentialOuterAlgorithmOptions, eps::Real)="comp_pen"

_compromise_cname(x) = _compromise_name(x) * "_fpath"
_compromise_fname(x) = _compromise_name(x) * ".jld2"

function run_compromise!(
    df::DF.DataFrame;
    settings = SettingsContainer(),
    colname = _compromise_cname(settings),
    kwargs...
)
    df_tmp = run_compromise(df; settings, colname, kwargs...)
    if colname in names(df)
        DF.select!(df, DF.Not(colname))
    end
    DF.leftjoin!(df, df_tmp; on=["problem_name", "constraint_index"])
    return df
end

function run_compromise(
    df::DF.DataFrame;
    settings = SettingsContainer(),
    colname = _compromise_cname(settings),
    thread_rows::Bool=false,
    disable_logging::Bool=false,
    kwargs...
)  
    
    @unpack num_halton, max_func_calls, max_points = settings
    if !_check_df_has_column(df, num_halton)
        error("DataFrame is missing column \"num_halton\" required by `settings`")
    end
    if !_check_df_has_column(df, max_func_calls)
        error("DataFrame is missing column \"max_func_calls\" required by `settings`")
    end
    if !_check_df_has_column(df, max_points)
        error("DataFrame is missing column \"max_points\" required by `settings`")
    end

    df_tmp = DF.DataFrame(
        "problem_name" => String[],
        "constraint_index" => Int[],
        colname => Union{Missing,String}[],
    )

    base_logger = if disable_logging 
        Logging.NullLogger()
    else
        Logging.current_logger()
    end

    disk_lock = ReentrantLock()
    
    return Logging.with_logger(base_logger) do
        run_compromise(
            Val(thread_rows), 
            df_tmp, df
            ;
            settings,
            disk_lock,
            kwargs...
        )
    end
end

_check_df_has_column(df, cname)=_check_df_has_column(df, string(cname))
_check_df_has_column(df, cname::Integer)=true
_check_df_has_column(df, cname::Nothing)=true
function _check_df_has_column(df, cname::AbstractString)
    return cname in names(df)
end

function run_compromise(
    thread_rows::Val{false}, 
    df_tmp::DF.DataFrame, df_lib::DF.DataFrame
    ;
    settings, 
    log_level = Info,
    logging_prefix = nothing,
    kwargs...
)
    @logmsg log_level "Non-Threaded Optimization."
    @progress for row in eachrow(df_lib)
        @unpack problem_name, constraint_index = row
        res_fpath = run_compromise_on_row(row, settings; log_level, kwargs...)
        push!(df_tmp, (problem_name, constraint_index, res_fpath))
    end
    GC.gc(false)
    return df_tmp
end

function run_compromise(
    thread_rows::Val{true},
    df_tmp::DF.DataFrame, df_lib::DF.DataFrame
    ;
    settings,
    log_level = Info,
    logging_prefix = nothing,
    kwargs...
)
    @logmsg log_level "Threaded Optimization."
    thread_logger = TransformerLogger(
        _make_log_interceptor(logging_prefix), 
        Logging.current_logger()
    )
    rlock = ReentrantLock()
    num_rows = size(df_lib, 1)  # for progress bar
    counter_rows = Threads.Atomic{Int}(0)        # for progress bar
    Logging.with_logger(thread_logger) do
        @withprogress begin
            Threads.@threads for row in eachrow(df_lib)
                @unpack problem_name, constraint_index = row
                res_fpath = run_compromise_on_row(row, settings; log_level, kwargs...)
                lock(rlock) do
                    push!(df_tmp, (problem_name, constraint_index, res_fpath))
                end
                
                Threads.atomic_add!(counter_rows, 1)
                @logprogress counter_rows[] / num_rows
            end
        end
    end
    GC.gc(false)
    return df_tmp
end
function _make_log_interceptor(::Nothing)
    _make_log_interceptor("")
end

function _make_log_interceptor(prefix)
    return Base.Fix1(_log_interceptor, prefix)
end

function _log_interceptor(prefix, log_args)
    @unpack level, message, _module, group, id, file, line, kwargs = log_args
    message = "$(prefix)$(Dates.format(Dates.now(), "HH:MM:SS.sss")) (THREAD $(Threads.threadid())) -- " * message
    return (; level, message, _module, group, id, file, line, kwargs)
end

function run_compromise_on_row(
    row::DF.DataFrameRow,
    user_settings :: Union{Nothing, SettingsContainer} = SettingsContainer() 
    ;
    results_name = _compromise_fname(user_settings),
    overwrite_results::Bool=false,
    check_settings::Bool=true,
    disk_lock,
    kwargs...
)
    
    @unpack problem_name, constraint_index, f90_fpath = row
   
    subdir = abspath(first(splitdir(f90_fpath)))
    res_fpath = joinpath(subdir, results_name)
    if overwrite_results 
        rm(res_fpath; force=true)
    end

    settings = _parse_settings_container(user_settings, row)
    _is_valid = lock(disk_lock) do
        has_valid_field(res_fpath, "settings", settings; check_data = check_settings)
    end
    if _is_valid
        return res_fpath
    end
    return maybe_run_compromise(
        row, res_fpath, settings; 
        disk_lock, kwargs...
    ) 
end

function _parse_settings_container(settings::Nothing, row)
    return nothing
end
function _parse_settings_container(settings::SettingsContainer, row)
    @unpack num_halton, max_func_calls, max_points = settings
    num_halton = _row_settings_val(row, num_halton)
    max_func_calls = _row_settings_val(row, max_func_calls)
    max_points = _row_settings_val(row, max_points)
    @unpack kernel, eps = settings
    algo_opts = deepcopy(settings.algo_opts)
    return SettingsContainer(;
        num_halton, max_func_calls, max_points,
        algo_opts, kernel, eps 
    )
end
_row_settings_val(row, cname)=_row_settings_val(row, string(cname))
_row_settings_val(row, cname::AbstractString)=row[cname]
_row_settings_val(row, cname::Integer)=cname
_row_settings_val(row, cname::Nothing)=cname

function maybe_run_compromise(row, res_fpath, settings::Nothing; kwargs...)
    GC.safepoint()
    return res_fpath
end

function maybe_run_compromise(
    row, res_fpath, parsed_settings::SettingsContainer; 
    disk_lock,
    log_level = Info,
    inner_log_level = nothing,
    kwargs...
)
    settings = _reset_log_level(parsed_settings, inner_log_level)
    ret_val = res_fpath
    @unpack problem_name, constraint_index, so_fpath = row
    @unpack num_halton = settings
    x0, lb, ub = BFH._row_dat(row, "x0", "lb", "ub")
    x0_matrix = hcat(x0, sample_halton(num_halton, lb, ub))
    res0 = values_for_matrix(so_fpath, x0_matrix)
    try 
        @logmsg log_level "[$(Threads.threadid())] Compromise on $(problem_name) / $(constraint_index)."
        x_opt, num_evals = run_compromise_on_matrix(so_fpath, x0_matrix, settings;)
        res_opt = values_for_matrix(so_fpath, x_opt)

        lock(disk_lock) do
            JLD2.jldsave(
                res_fpath
                ; 
                num_evals, settings = parsed_settings, 
                x0=res0.x, fx0=res0.fx, viol0=res0.viol, fx_eps0=res0.fx_eps, 
                res_opt...,
            )
        end
    catch err
        @error "Could not finish $(problem_name) / $(constraint_index):" exception=(err, catch_backtrace())
        ret_val = missing
    end
    return ret_val
end

_reset_log_level(settings::SettingsContainer, ::Nothing) = settings
function _reset_log_level(settings::SettingsContainer, log_level)
    @reset settings.algo_opts = _reset_log_level(settings.algo_opts, log_level)
    return settings
end
function _reset_log_level(algo_opts::SequentialOuterAlgorithmOptions, log_level)
    @reset algo_opts.log_level = log_level
    @reset algo_opts.inner_opts.log_level = log_level
    return algo_opts
end

function _reset_log_level(algo_opts::AlgorithmOptions, log_level)
    @reset algo_opts.log_level = log_level
    return algo_opts
end

function run_compromise_on_matrix(problem_so_path, x0_matrix, settings)
    @unpack algo_opts, kernel, max_func_calls, max_points, eps = settings
    return run_compromise_on_matrix(
        problem_so_path, x0_matrix, algo_opts, eps; 
        kernel, max_func_calls, max_points
    )
end

function run_compromise_on_matrix(
    problem_so_path, x0_matrix, algo_opts::SequentialOuterAlgorithmOptions, eps::Nothing;
    kernel, max_func_calls, max_points
)
    run_compromise_multi(
        problem_so_path, x0_matrix;
        algo_opts, kernel, max_func_calls, max_points
    )
end

function run_compromise_on_matrix(
    problem_so_path, x0_matrix, algo_opts::AlgorithmOptions, eps;
    kernel, max_func_calls, max_points
)
    return run_compromise_set(
        problem_so_path, x0_matrix;
        algo_opts, kernel, max_func_calls, max_points
    )
end

function run_compromise_on_matrix(
    problem_so_path, x0_matrix, algo_opts::SequentialOuterAlgorithmOptions, eps::Real;
    kernel, max_func_calls, max_points
)
   return run_compromise_penalized(
        problem_so_path, x0_matrix;
        algo_opts, kernel, max_func_calls, eps, max_points
    )
end

function run_compromise_multi(
    problem_so_path, x0_matrix;
    algo_opts=SequentialOuterAlgorithmOptions(), 
    kernel=CubicKernel(),
    max_func_calls=2_000,
    max_points=nothing
)
    return fortran_wrap(problem_so_path) do evaluator
        num_vars = evaluator.n_vars
        num_objfs = evaluator.n_objfs
        num_constrs = evaluator.n_constrs

        lb = copy(evaluator.lb)
        ub = copy(evaluator.ub)
        mop = MutableMOP(; 
            num_vars, lb, ub,
            reset_call_counters=false,
        )

        objfs_handle = inplace_objectives(evaluator)        
        objfs_cfg = surrogate_config(
            kernel, num_vars, num_objfs, max_points
        )
        add_objectives!(
            mop, objfs_handle, objfs_cfg;
            max_func_calls, dim_out=num_objfs, func_iip=true,
        )
        
        if num_constrs > 0
            constrs_handle = inplace_constraints(evaluator)
            constrs_cfg = surrogate_config(kernel, num_vars, num_constrs, max_points)
            add_nl_ineq_constraints!(
                mop, constrs_handle, constrs_cfg;
                max_func_calls, dim_out=num_constrs, func_iip=true,
            )
        end
        
        compromise_return_objects = optimize_with_algo(mop, algo_opts, x0_matrix)
       
        num_evals_objfs = evaluator.counter_objfs[]
        num_evals_constrs = evaluator.counter_constrs[]
        num_evals = max(num_evals_objfs, num_evals_constrs)

        X = parse_return_objects(compromise_return_objects, num_vars)
        return X, num_evals    
    end
end

function run_compromise_penalized(
    problem_so_path, x0_matrix;
    algo_opts=SequentialOuterAlgorithmOptions(), 
    kernel=CubicKernel(),
    max_func_calls=2_000,
    eps=0.1,
    max_points=nothing
)
    return fortran_wrap(problem_so_path) do evaluator
        num_vars = evaluator.n_vars
        num_objfs = evaluator.n_objfs

        lb = copy(evaluator.lb)
        ub = copy(evaluator.ub)
        mop = MutableMOP(; 
            num_vars, lb, ub,
            reset_call_counters=false,
        )
  
        evaluator.counter_objfs[] = 0
        evaluator.counter_constrs[] = 0
        
        mop_objfs_handle = l1_penalized_objectives(evaluator; eps)
        #mop_objfs_handle = CTS.lu_penalized_objectives(evaluator; eps)
        objfs_cfg = surrogate_config(kernel, num_vars, num_objfs, max_points)
        add_objectives!(
            mop, mop_objfs_handle, objfs_cfg;
            max_func_calls, dim_out=num_objfs, func_iip=true,
        )
        
        compromise_return_objects = optimize_with_algo(mop, algo_opts, x0_matrix)
        
        num_evals_objfs = evaluator.counter_objfs[]
        num_evals_constrs = evaluator.counter_constrs[]
        num_evals = max(num_evals_objfs, num_evals_constrs)

        X = parse_return_objects(compromise_return_objects, num_vars)
        return X, num_evals
    end
end

function run_compromise_set(
    problem_so_path, x0_matrix;
    algo_opts=Compromise.AlgorithmOptions(), 
    kernel=CubicKernel(),
    max_func_calls=2_000,
    max_points=nothing
)
    return fortran_wrap(problem_so_path) do evaluator
        num_vars = evaluator.n_vars
        num_objfs = evaluator.n_objfs
        num_constrs = evaluator.n_constrs        

        lb = copy(evaluator.lb)
        ub = copy(evaluator.ub)
        mop = MutableMOP(; 
            num_vars, lb, ub,
            reset_call_counters=false,
        )
        evaluator.counter_objfs[] = 0
        evaluator.counter_constrs[] = 0

        objfs_handle = inplace_objectives(evaluator)
        
        objfs_cfg = surrogate_config(kernel, num_vars, num_objfs, max_points)

        add_objectives!(
            mop, objfs_handle, objfs_cfg;
            max_func_calls, dim_out=num_objfs, func_iip=true,
        )
        if num_constrs > 0
            constrs_handle = inplace_constraints(evaluator)        
            constrs_cfg = surrogate_config(kernel, num_vars, num_constrs, max_points)
            add_nl_ineq_constraints!(
                mop, constrs_handle, constrs_cfg;
                max_func_calls, dim_out=num_constrs, func_iip=true,
            )
        end
        
        population = optimize_many(x0_matrix, mop; algo_opts)
        
        num_evals_objfs = evaluator.counter_objfs[]
        num_evals_constrs = evaluator.counter_constrs[]
        num_evals = max(num_evals_objfs, num_evals_constrs)

        return Compromise.cached_ξ(population), num_evals
    end
end

function parse_return_objects(compromise_return_objects, num_vars)
    num_opt = length(compromise_return_objects)
    X = zeros(num_vars, num_opt)
    keep_col = ones(Bool, num_opt)
    for (i, ret_obj) in enumerate(compromise_return_objects)
        val_cache = Compromise.opt_vals(ret_obj)
        if isnothing(val_cache)
            keep_col[i] = false
        else
            X[:, i] .= Compromise.cached_ξ(val_cache)
        end
    end
    X = X[:, keep_col]
    return X
end

function surrogate_config(kernel, num_vars, num_out, max_points=nothing)
    rbf_db = Compromise.RBFModels.init_rbf_database(num_vars, num_out, nothing, nothing, Float64)
    return Compromise.RBFModels.RBFConfig(; kernel, database=rbf_db, max_points)
end
#=
function sample_halton(row, settings)
    @unpack num_vars = row
    max_func_calls = settings_max_func_calls(settings, num_vars)
        
    n_halton = min(
        settings_num_halton(settings, num_vars),
        ceil(Int, max_func_calls / 3)
    )

    x0, lb, ub = BFH._row_dat(row, "x0", "lb", "ub")
    return hcat(x0, sample_halton(n_halton, lb, ub))
end
=#
function sample_halton(n_halton, lb, ub)
    num_vars = length(lb)
    w = ub .- lb
    return mapreduce(
        x -> lb .+ w .* x, 
        hcat, 
        HaltonPoint(num_vars; length=n_halton); 
        init=zeros(num_vars, 0)
    )
end

end