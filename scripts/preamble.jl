# # Common Definitons 

# This file contains code that is relevant to all (or most) of the scripts.

# ## Path Definitons
# Paths are defined in a separate file for lightweight inclusion.
include("path_definitions.jl")

# ## Dependencies
# We use `Pkg` to activate the environment in the directory containing this script.
# The environment is set up to provide all the necessary packages.
using Pkg
Pkg.activate(JULIA_PROJECT_ENV)

# We rely heavily on `DataFrames` for tabular data:
import DataFrames as DF
# Files are stored in CSV or JLD2 format.
import JLD2
# `MOBenchmarks` is a local package to help with running Fortran code and performing 
# analysis on optimization results.
import MOBenchmarks as MOB

# Our own optimizer is `Compromise`:
import Compromise as C
# We can also already load the Compromise benchmark extension:
CH = MOB.init_CompromiseHelpers()
# Currently, it “ships with” HiGHS for the linear or quadratic subproblems.
# If we wanted to, we could load `OSQP` and use it instead (as an extension).
# import OSQP

# The package `FiniteDiff` allows for easy finite-difference derivatives..
import FiniteDiff as FD

# For plotting, we use (GL)Makie and `LaTeXStrings`:
# using GLMakie
# import GLMakie: RGB, RGBA
using CairoMakie
import CairoMakie: RGB, RGBA
using LaTeXStrings

# With `@unpack` I can be lazy about named tuples:
import UnPack: @unpack
# Similarly, I use `@set` and `@reset` for immutable objects:
import Accessors: @set, @reset

# The standard `Dict` does not keept its order.
# I use an `OrderedDict` to store results.
import DataStructures: OrderedDict

# To disable logging, we need the corresponding package:
import Logging
import Logging: @logmsg, global_logger, Debug, Warn, Info

## To have progress logging available outside of VSCode:
import TerminalLoggers: TerminalLogger
if get(ENV, "TERM_PROGRAM", "") != "vscode"
    global_logger(TerminalLogger())
end

# A helper function to make sure that a base-path exists when saving files:
import MOBenchmarks.FortranHelpers: ensure_path

# ## Makie Theme
WONG_COLORS = Makie.wong_colors()
CUSTOM_THEME = merge(
    theme_latexfonts(),
    Theme(;
        fontsize = 10f0,        ## base fontsize; most attributes inherit this
        Axis = (;
            titlesize = 12f0,
            xlabelsize = 12f0,
            ylabelsize = 12f0,
        ),
        Lines = (;
            linewidth = 1.5f0,
        ),
        Scatter = (;
            markersize = 7f0,
            strokewidth = .75f0,
        ),
        ScatterLines = (;
            linewidth = 1.5f0,
            markersize = 7f0,
            strokewidth = .75f0,
        ),
        Legend = (;
            rowgap = 0f0,
            patchsize = (12f0, 12f0),
            padding = (3f0, 3f0, 3f0, 3f0)
        )
    ),
)
CUSTOM_THEME_SMALL = merge(
    CUSTOM_THEME,
    Theme(;
        Axis = (;
            xlabelpadding = 1f0,
            ylabelpadding = 1f0,
        ),
        Lines = (;
            linewidth = 1f0,
        ),
        Scatter = (;
            markersize = 2f0,
            strokewidth = .25f0,
        ),
        ScatterLines = (;
            linewidth = 1f0,
            markersize = 2f0,
            strokewidth = .25f0,
        ),
        Legend = (;
            padding = (1f0, 1f0, 1f0, 1f0)
        )
    ),
)

function dark_color(rgb::RGB)
    @unpack r, g, b = rgb
    return dark_color(r, g, b)
end
function dark_color(rgb::RGBA)
    @unpack r, g, b = rgb
    return RGBA(dark_color(r, g, b), rgb.alpha)
end
function dark_color(r, g, b)
    r, g, b = .5 .* (r, g, b)
    return RGB(r, g, b)
end