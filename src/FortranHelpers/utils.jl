
isnearly(a, b) = length(a) != length(b) ? false : isapprox(a, b)

function prepare_target_path(
    lib_path;
    log_level = Warn,
    indent = 0,
    overwrite=false, 
    return_if_exists=lib_path,
    return_if_fail=nothing,
    want_empty=true
)
    if ispath(lib_path)
        if isdir(lib_path)
            if !isempty(readdir(lib_path))
                if overwrite
                    rm(lib_path; recursive = true)
                else
                    if want_empty
                        @logmsg log_level "$(indent_str(indent))There is a non-empty directory at $(lib_path)."
                    end
                end
            end
        else
            if overwrite
                rm(lib_path)
            else
                if want_empty
                    @logmsg log_level "$(indent_str(indent))There is a file at $(lib_path)"
                end
            end
        end
    end
    _lib_path = return_if_exists
    if !ispath(lib_path)
        _lib_path = return_if_fail
        try
            mkpath(lib_path)
            _lib_path = lib_path
        catch err
            @error "Could not make path `lib_path`." exception = (err, catch_backtrace())
        end
    end
    return _lib_path
end

function data_from_shared_lib(so_path)
    x0, lb, ub, n_vars, n_objfs, n_constrs = fortran_wrap(so_path) do evaluator
        return (
            copy(evaluator.x0), 
            copy(evaluator.lb), 
            copy(evaluator.ub),
            evaluator.n_vars, 
            evaluator.n_objfs, 
            evaluator.n_constrs
        )
    end
    return x0, lb, ub, n_vars, n_objfs, n_constrs
end

function _row_dat(row, fields...)
    return JLD2.load(row["dat_fpath"], fields...)
end

function _row_feas(row, fields...)
    return JLD2.load(row["feas_fpath"], fields...)
end

augment_dat!(df) = _augment!(df, "dat_fpath", "x0", "lb", "ub")
augment_feas!(df) = _augment!(df, "feas_fpath", "xfeas", "x0_feasible", "theta", "nlopt_ret")

augment_dat(df) = _augment(df, "dat_fpath", "x0", "lb", "ub")
augment_feas(df) = _augment(df, "feas_fpath", "xfeas", "x0_feasible", "theta", "nlopt_ret")

function _augment(df::DF.DataFrame, colname, fields...)
    return _augment_with_func(DF.select, df, colname, fields...)
end
function _augment!(df::DF.DataFrame, colname, fields...)
    return _augment_with_func(DF.select!, df, colname, fields...)
end

function _augment_with_func(selector, df, colname, fields...)
    if colname in names(df)
        new_cols = collect(fields)
        loader = dat_fpath -> JLD2.load(dat_fpath, fields...)
        _df = selector(
            df, 
            :,
            colname => DF.ByRow( loader ) => new_cols
        )
    else
        _df = df
    end
    return _df
end

_intify(i::Real)=Int(i)
_intify(i::String)=parse(Int, i)

indent_str(::Nothing) = ""
indent_str(i::Integer) = lpad("", i)