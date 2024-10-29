# Setup Script

# This script sets up the Julia environment according to `JULIA_PROJECT_ENV` as in 
# `path_definitions`:
include("path_definitions.jl")

# Delete old environment files for a fresh install:
rm(JULIA_PROJECT_ENV; force=true, recursive=true)

# Activate new environment:
using Pkg;
Pkg.activate(JULIA_PROJECT_ENV);
# Add `Compromise` and Benchmark helpers.
# If `Pkg.develop` is used, we can benefit from `Revise.jl` and code reloading more easily.
Pkg.develop(; url="https://github.com/manuelbb-upb/Compromise.jl.git")
Pkg.develop(; path=BENCHMARK_PACKAGE_PATH)
# Packages for plotting:
Pkg.add([
    "GLMakie",
    "LaTeXStrings",
])
# Other helpers:
Pkg.add([
    "Accessors",
    "UnPack",
    "DataFrames",
    "DataStructures",
    "JLD2",
    "TerminalLoggers"
])
# Other optimization-related packages:
Pkg.add("NLopt")
Pkg.add("FiniteDiff")
Pkg.add("HaltonSequences")
# Pkg.add("OSQP")   # see explanation in `preamble.jl`