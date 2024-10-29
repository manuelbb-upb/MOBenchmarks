# Introduction 
This repository contains code to compare the derivative-free, multi-objective 
solvers `DFMO` and `Compromise.jl` on a set of nonlinearly constrained test problems 
derived from `TESTMO`.
Both `DFMO` and `TESTMO` are written in Fortran.
`Compromise` is written in Julia and provides various algorithms
for multi-objective optimization.
It turns out that benchmarking optimizers written in different languages on a common
set of test problems is a bit messy.
The messiness is readily reflected in code.
We have tried to keep a semblance of order by abstracting the most 
technical aspects, e.g., compilation of external libraries, 
and by separating chunks of codes into modules.
Everything is encapsulated in a big module `MOBenchmarks`. 
It has a `FortranHelpers` submodule.

As `Compromise.jl` is not registered in the general registry (yet),
it is not a strict dependency of `MOBenchmarks`, but it enables an **extension** 
if present.

# Tabular Raw Data

We take (a subset) of problems from `TESTMO` and augment them by adding nonlinear constraints.
A test problem is then defined by a tuple of problem name and constraint index.
A constraint index of `0` refers to the original, unmodified problem.
Other possible constraint indices range from `1` to `6`.

Most of the time we use `DataFrames.DataFrame` objects for “organizing” and filtering
the test problems.
Some operations (feasibility testing and optimization) require much time.
That is why we store results in files.
For simple tabular data, we sometimes use the `CSV` file format.
Otherwise, `JLD2` is used to store Julia objects.
This way, we can also store configuration objects of custom types,
and use those to check if an operation has to be redone.

# Scripts
How to use the `scripts`.

* First, make sure that `path_definitions.jl` suits your needs.
* Run `setup_julia.jl` to install the required Julia packages.
* The problem library is instantiated with `performance_run_dfmo.jl`.
* Afterwards, `performance_run_compromise.jl` can be used.
* Finally, `performance_postprocess.jl` generates result plots.
* The files `example_mw3.jl` and `example_two_paraboloids.jl` are separate from the  
  performance analysis.

Beware, that subsequent execution of the scripts from the same REPL might clutter the 
global workspace and lead to unexpected results.
A possible workaround is to modify the scripts to use modules.
Or `include` within a `let` scope.

# Dev Environment

In an attempt to declutter my global environment, I have tried out the devcontainer feature of VSCode.
For various reasons, it did not work out for me, but I have kept the `.devcontainer` folder.
Instead, I now rely on `flake.nix` to have a more isolated and reproducible environment.