@enum DFMO_DIRTYPE :: UInt8 begin
    DFMO_HALTON_DIRS = 1
    DFMO_SOBOL_DIRS = 2
end

"""
    DFMOSettings(; kwargs...)

Return a mutable settings object for compilation of a problem for DFMO.

# Keyword arguments

* `xopt0 :: Union{Vector{Float64}, Missing} = missing`: starting vector for optimization.  
  If `missing`, then the `startp` subroutine is called.
* `alfa_stop :: Float64 = 1e-9`: minimum step length.
* `nf_max :: Integer = 2000`: maximum number of function evaluations.
* `iprint :: Integer = 0`: printing level. A value of `0` disables console output.
* `hschoice :: DFMO_DIRTYPE = DFMO_SOBOL_DIRS`: which directions to use.
* `dir_dense :: Bool = true`: whether to use dense direction or not.
* `dir_coord :: Bool = true`: whether to use coordinate directions or not.
"""
Base.@kwdef mutable struct DFMOSettings
    xopt0 :: Union{Vector{Float64}, Missing} = missing
    alfa_stop :: Float64 = 1e-9   # tolerance for step_length termination. 
                            # DFMO will terminate as soon as all of the step_length 
                            # fall below alfa_stop
    nf_max :: Integer = 2_000       # maximum number of allowed function evaluations
    iprint :: Integer = 0           # printing level. 0 - no console output, >0 different levels of printing
    hschoice :: DFMO_DIRTYPE = DFMO_SOBOL_DIRS    # which type of dense direction is used. 1 - HALTON-type, 2 - SOBOL-type
    dir_dense :: Bool = true        # whether to use the dense direction or not
    dir_coord :: Bool = true        # whether to use the coordinate directions or not
end

@batteries DFMOSettings selfconstructor=false

const DEFAULT_TMP_DIR = tempname()

const DFMO_FILES = (
    "modules_DFMO",
    "main",
    "DFMO",
    "subroutines_DFMO",
    "halton",
    "sobol",
    "qsortd",
    "problem",
)

function compile_dfmo(
    problem_src_path, 
    output_file_path, 
    dfmo_src_dir,
    dfmo_settings :: DFMOSettings = DFMOSettings();
    build_path= DEFAULT_TMP_DIR,
    remove_build_dir::Bool= false,
    force_recompile::Bool=true,
    log_level=Info,
    logmsg_build=nothing
)
    @unpack xopt0, alfa_stop, iprint, nf_max, hschoice, dir_dense, dir_coord = dfmo_settings
    hschoice = Int(hschoice)
    return compile_dfmo(
        problem_src_path, output_file_path, dfmo_src_dir,
        xopt0, alfa_stop, nf_max, iprint, hschoice, dir_dense, dir_coord;
        build_path, remove_build_dir, force_recompile,
        log_level, logmsg_build
    )
end

function compile_dfmo(
    problem_src_path, 
    output_file_path, 
    dfmo_src_dir,
    xopt0::Union{Vector{Float64}, Missing},
    alfa_stop::Float64,
    nf_max::Integer,
    iprint::Integer,
    hschoice::Integer,
    dir_dense::Bool,
    dir_coord::Bool;
    build_path=DEFAULT_TMP_DIR,
    remove_build_dir::Bool=false,
    force_recompile::Bool=true,
    log_level=Info,
    logmsg_build=nothing
)
    output_file_path = abspath(output_file_path)

    if !isfile(output_file_path) || force_recompile
        if isa(logmsg_build, AbstractString)
            @logmsg log_level logmsg_build
        end
        cp_trgt_path = joinpath(build_path, "DFMO")
     
        mkpath(cp_trgt_path)
        for fn in readdir(dfmo_src_dir)
            if endswith(fn, ".f90") || endswith(fn, ".inc") || fn == "makefile"
                cp(joinpath(dfmo_src_dir, fn), joinpath(cp_trgt_path, fn); force=true)
            end
        end

        dfmo_hashes_file_path = joinpath(cp_trgt_path, "src_hashes.jldat")
        dfmo_hashes = isfile(dfmo_hashes_file_path) ? deserialize(dfmo_hashes_file_path) : Dict{String, String}()
        
        trgt_f90_path = joinpath(cp_trgt_path, "problem.f90")
        cp(problem_src_path, trgt_f90_path; force=true)
        if !ismissing(xopt0)
            change_startp_in_problem_f90!(trgt_f90_path, xopt0)
        end

        main_path = joinpath(cp_trgt_path, "main.f90")
        open(main_path, "r+") do file
            main_content = read(file, String)
            main_content = replace(
                main_content,
                r"(\N*?alfa_stop\s*?=\s*)([\w\.\-]*)"m => SubstitutionString("\\g<1>" * @sprintf("%.15f", alfa_stop)),
                r"(\N*?nf_max\s*?=\s*)([\w\.\-]*)"m => SubstitutionString("\\g<1>" * "$(nf_max)"),
                r"(\N*?iprint\s*?=\s*)([\w\.\-]*)"m => SubstitutionString("\\g<1>" * "$(iprint)"),
                r"(\N*?hschoice\s*?=\s*)([\w\.\-]*)"m => SubstitutionString("\\g<1>" * "$(hschoice)"),
                r"(\N*?dir_dense\s*?=\s*)([\w\.\-]*)"m => SubstitutionString("\\g<1>" * (dir_dense ? ".true." : ".false.")),
                r"(\N*?dir_coord\s*?=\s*)([\w\.\-]*)"m => SubstitutionString("\\g<1>" * (dir_coord ? ".true." : ".false.")),
            )
            truncate(file, 0)
            write(file, main_content)
        end

        current_dir = pwd()
        try
            cd(cp_trgt_path)
            #run(`make`)
            for fn in DFMO_FILES
                fn_old_hash = get(dfmo_hashes, fn, "")
                fn_src_path = joinpath(cp_trgt_path, "$(fn).f90")
                fn_trgt_path = joinpath(cp_trgt_path, "$(fn).o")
                fn_hash = hash_file_contents(fn_src_path)
                if fn_hash != fn_old_hash
                    dfmo_hashes[fn] = fn_hash
                    run(`gfortran -O3 -c $(fn_src_path) -o $(fn_trgt_path)`)
                end
            end
            run(`gfortran -o multiobj modules_DFMO.o main.o DFMO.o subroutines_DFMO.o halton.o sobol.o qsortd.o problem.o`)
            cp("multiobj", output_file_path; force=true)
            serialize(dfmo_hashes_file_path, dfmo_hashes)
        catch err
            rethrow(err)
        finally
            cd(current_dir)
        end

        if remove_build_dir
            rm(cp_trgt_path; recursive=true)
        end
        
    end
    return nothing
end