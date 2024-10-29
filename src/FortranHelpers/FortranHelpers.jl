module FortranHelpers

## imports from standard library
import Printf: @sprintf         # for simple manipulation (substitution) in Fortran source code
import Libdl                    # working with shared library objects (compiled Fortran)
import Serialization: serialize, deserialize    # serialize and store hashes of source code files to cache dependencies in Fortran compilation
import SHA: sha2_256            # hash file contents
import Logging: @logmsg, Info, Warn, LogLevel

## dependencies for downloading DFMO and TESTMO
import Downloads as DL
import Tar
import CodecZlib: GzipDecompressorStream

## helpers for structs and named tuples
import Parameters: @with_kw
import UnPack: @unpack
import StructHelpers: @batteries

## tabular data management
import DataFrames as DF
import DataFrames: DataFrame, metadata

## storing files
import JLD2
import CSV

## nonlinear optimization
import NLopt

## progress bar in compatible loggers and bottom of vscode terminal
import ProgressLogging: @progress

include("constants.jl")

include("download_fortran.jl")
include("problem_lib.jl")
include("feasibility_checks.jl")
include("run_dfmo.jl")
#include("metrics.jl")

function hash_file_contents(file_path)
    return open(ensure_path(file_path)) do f
        bytes2hex(sha2_256(f))
    end
end

function ensure_path(p)
    base, ext = splitext(p)
    dir = isempty(ext) ? base : first(splitdir(p))
    if !isdir(dir)
        mkpath(dir)
    end
    return p
end

include("ineq_constraints.jl")
include("fortran_code_manipulation.jl")
include("dfmo_compilation.jl")
include("problem_evaluator.jl")
include("utils.jl")
export DFMOSettings

end # module CompromiseTestSuite
