# With this script we download the external dependencies `DFMO` and `TESTMO` to
# then generate the augmented test problem library at `LIB_PATH`.
# We also check the problems for feasibility and run `DFMO`.

# Import global definitions:
include("preamble.jl")

# Download external packages from GitHub:
MOB.download_dfmo(DFMO_PATH; overwrite=false)
MOB.download_testmo(TESTMO_PATH; overwrite=false)

# Scan `TESTMO` for `.f90` Fortran files and build a basic dataframe:
testmo_src_df = MOB.scan_testmo_dir(
    TESTMO_PATH; 
    log_level = Info
)

# From this, build the augmented problem library at `LIB_PATH` and return a dataframe 
# with metadata:
lib_base_df = MOB.initialize_problem_lib(
    LIB_PATH, testmo_src_df;
    ## keep the folder at `LIB_PATH` if it exists already:
    overwrite_problem_lib = false,
    ## keep the individual problem files if they exist
    overwrite_problem_files = false,
    ## keep shared library files if they exist
    force_recompile = false,
    log_level = Info
)

# Discard problems that have fewer than 3 variables …
lib_base_df = lib_base_df[lib_base_df.num_vars .>= 3, :]
# … and then modify the Fortran code to have constraints:
lib_df = MOB.make_problems_with_constraints(lib_base_df)

# We use `NLopt` to check for feasiblity of the problems.
# Behind the scenes, we group problems by their variable bounds 
# and nonlinear constraints to speed this test up.
MOB.check_feasibility!(
    lib_df; 
    ## include only feasible test problems in the dataframe
    only_feasible=true,
    # use standard settings for the tests
    feasibility_test_settings = MOB.FortranHelpers.FeasibilityTestSettings(),
    ## keep feasibility information files if they already exist
    overwrite_feasibility_info_files = false,
    ### if there already is information for a feasibility proxy problem, keep it
    overwrite_feasibility_subdir = false,
    overwrite_feasibility_subdir_files = false,
    ## column name that is added to the dataframe to point to the feasibility info file
    cnames = ["feas_fpath",],
    log_level = Info
)

# Now we run DFMO on the feasible problems.
# This will ad a column to `lib_df` pointing to the result files (which have JLD2
# format and some standardized fields):
MOB.run_dfmo!(
    lib_df, DFMO_PATH;
    colname = COL_NAME_DFMO,
    ## use standard settings for the optimization runs
    settings = MOB.FortranHelpers.DFMOSettings(),
    ## validate settings; if they are different from what is stored, re-run the experiments
    check_settings = true,
    ## generally, try to keep results and binary
    overwrite_results = false,
    overwrite_exec = false,
    log_level = Info,
)

# Finally, we store the Dataframe as a CSV.
# We only store columns with non-float values, so that there are no rounding errors.
# Such meta-data can be retrieved from files that are referenced in the columns
# by means of helper functions.
problem_lib_csv_fpath = MOB.store!(lib_df; outname=PROBLEM_LIB_CSV_FNAME) |> abspath
@assert problem_lib_csv_fpath == PROBLEM_LIB_CSV_FPATH
@logmsg Info "Stored library CSV at \n`$(PROBLEM_LIB_CSV_FPATH)`"