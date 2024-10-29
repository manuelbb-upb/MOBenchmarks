module MOBenchmarks

const __revise_mode__ = :eval

import ProgressLogging: @progress
import CSV
import DataFrames as DF
import JLD2
import UnPack: @unpack
import Reexport: @reexport
include("FortranHelpers/FortranHelpers.jl")

import .FortranHelpers: hash_file_contents
import .FortranHelpers: 
    download_dfmo, download_testmo, 
    ## problem library instantiation
    scan_testmo_dir, initialize_problem_lib, make_problems_with_constraints, 
    ### feasibility testing
    check_feasibility, check_feasibility!,
    ### utils 
    store!, readlib,
    augment_dat!, augment_feas!,
    ## DFMO runs
    run_dfmo, run_dfmo!

import Metaheuristics.PerformanceIndicators: hypervolume
import Statistics: median, quantile
include("metrics.jl")

function init_CompromiseHelpers()
    ext = Base.get_extension(@__MODULE__, :CompromiseHelpersExt)
    if isnothing(ext)
        @warn "Could not load extension. Is Compromise available?"
    end
    return ext
end

function first_row(df, problem_name, constraint_index)
    _df = df[ (df.problem_name .== problem_name) .& (df.constraint_index .== constraint_index), :]
    if isempty(_df)
        return nothing
    else
        return first(_df)
    end
end

export FortranHelpers
export init_CompromiseHelpers

end # module MOBenchmarks
