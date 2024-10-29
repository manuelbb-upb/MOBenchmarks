# This file contains functions for setting up the test problem library.

# # First Stage -- Parse Fortran Problem Source
@doc """
   scan_testmo_dir(testmo_path; log_level=Info, indent=nothing)

Scan directory `testmo_path` (or subdir `PROBLEMS`) 
for files ending in ".f90" and return a `DataFrame` with columns 
"problem_name", "src_fpath".
If the path is not valid, return `nothing`.

The dataframe has a metadata field "testmo_path" of style `:note` containing
the full path to the `PROBLEMS` folder."""
function scan_testmo_dir(testmo_path; log_level=Info, indent=nothing)
    @logmsg log_level "$(indent_str(indent))Scanning TESTMO path '$(testmo_path)'..."
    !isa(testmo_path, AbstractString) && return nothing
    testmo_path = abspath(testmo_path)
    if !(last(splitdir(testmo_path)) == "PROBLEMS")
        testmo_path = joinpath(testmo_path, "PROBLEMS")
    end
    !isdir(testmo_path) && return nothing
    df = _empty_df_scan()
    for fn in readdir(testmo_path)
        !endswith(fn, ".f90") && continue
        problem_name = first(splitext(fn))
        src_fpath = joinpath(testmo_path, fn)
        push!(df, (; problem_name = problem_name, src_fpath); )
    end
    DF.metadata!(df, "testmo_path", testmo_path; style=:note)
    return df
end

function _empty_df_scan()
    return DataFrame(
        "problem_name" => String[],
        "src_fpath" => String[],
    )
end

@doc """
    initialize_problem_lib(
        problem_lib_path,
        testmo_src_df,
        ;
        overwrite = false,
        force_recompile = overwrite,
        problem_f90_filename = PROBLEM_F90_FILENAME,
        problem_so_filename = PROBLEM_SO_FILENAME,
        problem_dat_filename = PROBLEM_DAT_FILENAME,
        log_level = Info,
        indent = 0,
        kwargs...
    )

Initialize the test problem library at `problem_lib_path` based on
information in `testmo_src_df`.
The dataframe `testmo_src_df` is expected to have columns `problem_name` and `src_fpath`.
For each row, we copy problem source file and create a shared library object.
The shared library object is loaded to create a data file.
The data file provides problem information such as boundaries, number of variables, etc.

Return a dataframe with columns 
```
[
    "problem_name", "num_vars", "num_objectives", "x0", "lb", "ub",
    "f90_fpath", "so_fpath", "dat_fpath"
]
```

The keyword arguments specify whether to delete the library path or if shared 
library objects should be overcritten.
They also define the naming of the files that are generated."""
function initialize_problem_lib(
    problem_lib_path,
    testmo_src_df,
    ;
    overwrite_problem_lib = false,
    overwrite_problem_files = overwrite_problem_lib,
    force_recompile = overwrite_problem_files,
    problem_f90_filename = PROBLEM_F90_FILENAME,
    log_level = Info,
    indent = 0,
    kwargs...
)
    problem_lib_path = prepare_target_path(
        problem_lib_path; 
        overwrite=overwrite_problem_lib, indent = indent+1
    )
    isnothing(problem_lib_path) && return problem_lib_path
    problem_lib_path = abspath(problem_lib_path)
    
    df = _empty_df_init_problem_lib()
    istr = indent_str(indent)
    @progress for row in eachrow(testmo_src_df)
        @unpack src_fpath, problem_name = row
        problem_dir = joinpath(problem_lib_path, problem_name)
        if !isdir(problem_dir)
            mkpath(problem_dir)
        end
        
        problem_f90_fpath = joinpath(problem_dir, problem_f90_filename)
        ischanged_f90_fpath = force_recompile
        if overwrite_problem_files || !isfile(problem_f90_fpath)
            @logmsg log_level "$(istr)Problem '$(problem_name)', copying src to\n$(istr)`$(problem_f90_fpath)`"
            cp(src_fpath, problem_f90_fpath; force=true)
            ischanged_f90_fpath = true
        end
        
        problem_attrs = _make_shared_lib_and_data_file(
            problem_dir, problem_f90_fpath; 
            log_level, indent, force_recompile=ischanged_f90_fpath, problem_f90_filename, kwargs...
        )
        if isnothing(problem_attrs)
            @logmsg log_level "$(istr)!Skipping problem '$(problem_name)', could not read attributes."
            continue
        end

        f90_fpath = problem_f90_fpath
        push!(df, (; problem_name, f90_fpath, problem_attrs... ); cols=:intersect)
    end

    DF.metadata!(df, "problem_lib_path", problem_lib_path; style=:note)
    DF.metadata!(df, "testmo_path", DF.metadata(testmo_src_df, "testmo_path", nothing); style=:note)
    
    return df
end

function _empty_df_init_problem_lib()
    return DataFrame(
        "problem_name" => String[],
        "num_vars" => Int[],
        "num_objectives" => Int[],
        #"x0" => Vector{Float64}[],
        #"lb" => Vector{Float64}[],
        #"ub" => Vector{Float64}[],
        "f90_fpath" => String[],
        "so_fpath" => String[],
        "dat_fpath" => String[],
    )
end

@doc """
    _make_shared_lib_and_data_file(
        problem_dir,
        problem_f90_fpath=nothing;
        <keywoard arguments>
    )

Helper function to create shared library object and data file.
If successful, a NamedTuple with problem attributes is returned.
Otherwise, `nothing` is returned."""
function _make_shared_lib_and_data_file(
    problem_dir,
    problem_f90_fpath=nothing; 
    force_recompile=false,
    problem_f90_filename = PROBLEM_F90_FILENAME,
    problem_so_filename = PROBLEM_SO_FILENAME,
    problem_dat_filename = PROBLEM_DAT_FILENAME,
    log_level=Info, 
    indent=0,
)
    ## make sure that the Fortran source file is a path
    if isnothing(problem_f90_fpath)
        problem_f90_fpath = joinpath(problem_dir, problem_f90_filename)
    end
    if isnothing(problem_f90_fpath)
        @logmsg log_level "$(indent_str(indent))Cannot find source file `$(problem_f90_fpath)`."
        return nothing
    end
    ## create shared library object
    problem_so_fpath = joinpath(problem_dir, problem_so_filename)
    isnew_problem_so_fpath = false
    if force_recompile || !isfile(problem_so_fpath)
        isnew_problem_so_fpath = true
        @logmsg log_level "$(indent_str(indent))Compiling to $(problem_so_fpath)."
        run(`gfortran $(problem_f90_fpath) -fno-underscoring -shared -o $(problem_so_fpath)`)
    end
    if !isfile(problem_so_fpath)
        @logmsg log_level "$(indent_str(indent))Compiled shared library object not found at `$(problem_so_fpath)`"
    end
    ## make and read data file
    problem_dat_fpath = joinpath(problem_dir, problem_dat_filename)
    _problem_attrs = _shared_lib_to_data_file(
        problem_dat_fpath, problem_so_fpath; 
        force_recreate=isnew_problem_so_fpath,
        log_level, indent
    )
    return (; so_fpath=problem_so_fpath, dat_fpath=problem_dat_fpath, _problem_attrs...)
end

function _shared_lib_to_data_file(
    data_fpath, so_fpath; 
    force_recreate=false,
    log_level=Info,
    indent=0
)
    if force_recreate
        rm(data_fpath; force=true)
    end
    if !isfile(data_fpath)
        x0, lb, ub, num_vars, num_objectives, num_constraints = data_from_shared_lib(so_fpath)
        JLD2.jldsave(data_fpath; x0, lb, ub, num_vars, num_objectives, num_constraints)
    else
        x0, lb, ub, num_vars, num_objectives, num_constraints = JLD2.load(
            data_fpath, "x0", "lb", "ub", "num_vars", "num_objectives", "num_constraints")
    end

    return (; x0, lb, ub, num_vars, num_objectives, num_constraints)
end

# # Second Stage -- Adding Constraints
function make_problems_with_constraints(
    problem_lib_base_df,
    problem_lib_path=nothing;
    overwrite_subdirs=false,
    force_recreate=overwrite_subdirs,
    log_level=Info,
    indent=0
)
    istr = indent_str(indent)
    if isnothing(problem_lib_path)
        problem_lib_path = DF.metadata(problem_lib_base_df, "problem_lib_path", nothing)
    end
    if isnothing(problem_lib_path)
        @logmsg log_level "$(istr)Problem library path `problem_lib_path` not valid."
        return nothing
    end

    @logmsg log_level "$(istr)Creating problems with constraints."
    df = _empty_df_lib_with_constraints()
    @progress for row = eachrow(problem_lib_base_df)
        @unpack problem_name, dat_fpath, so_fpath = row
        base_f90_fpath = row["f90_fpath"]
        
        problem_base_dir, f90_fname = splitdir(base_f90_fpath)
        dat_fname = last(splitdir(dat_fpath))
        so_fname = last(splitdir(so_fpath))
        num_vars, num_objectives = JLD2.load(dat_fpath, "num_vars", "num_objectives")

        for constraint_index in 0:6
            subdir = joinpath(problem_base_dir, string(constraint_index))
            subdir = prepare_target_path(
                subdir; want_empty=false, overwrite=overwrite_subdirs, indent=indent+1)
            isnothing(subdir) && continue
            subdir = abspath(subdir)
            f90_fpath = joinpath(subdir, f90_fname)
            isnew_f90_fpath = force_recreate
            if !isfile(f90_fpath) || force_recreate
                isnew_f90_fpath = true
                cp(base_f90_fpath, f90_fpath; force=true)
                change_dims_in_problem_f90!(
                    f90_fpath, num_vars, num_objectives, Val(constraint_index))
            end

            problem_attrs = _make_shared_lib_and_data_file(
                subdir, f90_fpath;
                problem_so_filename = so_fname,
                problem_dat_filename = dat_fname,
                force_recompile = isnew_f90_fpath,
                log_level, indent
            )
            isnothing(problem_attrs) && continue
            push!(
                df, 
                (;problem_name, constraint_index, f90_fpath, problem_attrs...); 
                cols=:intersect
            )
        end
    end

    for metakey in ("problem_lib_path", "testmo_path")
        DF.metadata!(
            df, 
            metakey,
            DF.metadata(problem_lib_base_df, metakey, nothing); 
            style=:note
        )
    end
    
    return df    
end

function _empty_df_lib_with_constraints()
    return DataFrame(
        "problem_name" => String[],
        "constraint_index" => Int[],
        "num_vars" => Int[],
        "num_objectives" => Int[],
        "num_constraints" => Int[],
        #"x0" => Vector{Float64}[],
        #"lb" => Vector{Float64}[],
        #"ub" => Vector{Float64}[],
        "f90_fpath" => String[],
        "so_fpath" => String[],
        "dat_fpath" => String[],
    )
end

function store!(
    lib_df :: DF.DataFrame;
    outdir = nothing,
    outname = "problem_lib.csv",
    relative_paths = true
)
    if isnothing(outdir)
        outdir = DF.metadata(lib_df, "problem_lib_path", missing)
        return store!(lib_df; outdir, outname)
    end
    if ismissing(outdir)
        error("`outdir` invalid.")
    end

    if !ispath(outdir)
        @warn "`$(outdir)` is not a path yet."
        mkpath(outdir)
    end
    out_fpath = joinpath(outdir, outname)

    if relative_paths
        lib_df = copy(lib_df)
        _lib_df = DF.dropmissing(lib_df, :; view=true)
        for cname in names(_lib_df)
            if endswith(cname, "path")
                _lib_df[!, cname] .= relpath.(_lib_df[!, cname], Ref(outdir))
            end
        end
    end
    CSV.write(out_fpath, lib_df)
    return out_fpath
end

function readlib(
    csv_fpath;
    problem_lib_path = nothing,
    absolute_paths = true,
)
    lib_df = DF.DataFrame(CSV.File(csv_fpath))
    if isnothing(problem_lib_path)
        problem_lib_path = abspath(first(splitdir(csv_fpath)))
    end
    if absolute_paths
        for cname in names(lib_df)
            if endswith(cname, "path")
                lib_df[!, cname] .= map(
                    v -> ismissing(v) ? v : joinpath(problem_lib_path, v), lib_df[!, cname])
            end
        end
    end

    DF.metadata!(lib_df, "problem_lib_path", problem_lib_path; style=:note)
    return lib_df
end