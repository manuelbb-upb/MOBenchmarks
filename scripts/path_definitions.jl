# # Paths
# Path to generate files and folders in:
const BASE_DIR = joinpath(@__DIR__, "..", "results") |> abspath

# Julia project environment path:
const JULIA_PROJECT_ENV = joinpath(BASE_DIR, "JULIA_ENV") |> abspath

# Path to benchmark helper package:
const BENCHMARK_PACKAGE_PATH = joinpath(BASE_DIR, "..") |> abspath

# Path to `DFMO` solver source
const DFMO_PATH = joinpath(BASE_DIR, "DFMO") |> abspath
# Path to `TESTMO` test problem library:
const TESTMO_PATH = joinpath(BASE_DIR, "TESTMO") |> abspath
# Path to augmented test problem library.
# Results are also stored here:
const LIB_PATH = joinpath(BASE_DIR, "PROBLEMLIB") |> abspath

# Name of and path to base library CSV file:
const PROBLEM_LIB_CSV_FNAME = "problem_lib.csv"
const PROBLEM_LIB_CSV_FPATH = joinpath(LIB_PATH, PROBLEM_LIB_CSV_FNAME) |> abspath

# Names for tables containing Compromise results:
const LIB_SUFFIX = "_cubic_lp"
const LIB_CONSTR_CSV_NAME = "lib_constr$(LIB_SUFFIX).csv"
const LIB_UNCONSTR_CSV_NAME = "lib_unconstr$(LIB_SUFFIX).csv"
const LIB_CONSTR_CSV_FPATH = joinpath(LIB_PATH, LIB_CONSTR_CSV_NAME) |> abspath
const LIB_UNCONSTR_CSV_FPATH = joinpath(LIB_PATH, LIB_UNCONSTR_CSV_NAME) |> abspath

PLOTS_PATH = joinpath(BASE_DIR, "PLOTS") |> abspath
plots_path = PLOTS_PATH

# # Column Names
const COL_NAME_DFMO = "dfmo_fpath"
const COL_NAME_COMPROMISE_SET = "cset_fpath"
const COL_NAME_COMPROMISE_SEQ = "cseq_fpath"
const COL_NAME_COMPROMISE_PEN = "cpen_fpath"