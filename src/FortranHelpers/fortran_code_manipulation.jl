import Printf: @sprintf

const SET_DIMS_REGEX = r"subroutine setdim((.|\n)*)end subroutine setdim"
const STARTP_REGEX = r"subroutine startp((.|\n)*)end subroutine setdim"

function setdims_str(num_vars, num_nl_ineq_constrs, num_objfs)
    return """
    subroutine setdim(n, m, q)
        implicit none
        integer	:: n, m, q

        n = $(num_vars)
        m = $(num_nl_ineq_constrs)
        q = $(num_objfs)

        return
    end subroutine setdim"""
end

function startp_str(x0)
    x0_str = join((@sprintf("%.17f", xi) for xi in x0), ", ")
    return """
    subroutine startp(n,x)
        implicit none
        integer	:: n
        real*8	:: x(n)

        x = (/$(x0_str) /)
        
        return
    end subroutine startp"""
end

function change_startp_in_problem_f90!(trgt_file, x0)
    open(ensure_path(trgt_file), "r+") do f90_file
        prob_str = read(f90_file, String)
        ## change `startp` to return custom starting point vector
        new_prob_str = replace(
            prob_str, 
            STARTP_REGEX => startp_str(x0)
        )
        truncate(f90_file, 0)
        write(f90_file, new_prob_str)
    end
end

function change_dims_in_problem_f90!(
    trgt_file, n_vars, n_objfs, v::Val{constraint_index}
) where constraint_index
    n_constrs = num_vars_to_num_constraints(v, n_vars)
    open(ensure_path(trgt_file), "r+") do f90_file
        prob_str = read(f90_file, String)
        ## first, change `setdim` to return correct number of constraints
        new_prob_str = replace(
            prob_str, 
            SET_DIMS_REGEX => setdims_str(n_vars, n_constrs, n_objfs)
        )
        ## additionally, write the constraint function definition to the file:
        if !isnothing(constraint_index)
            new_prob_str *= "\n" * fortran_constraint_string(v)
        end
        truncate(f90_file, 0)
        write(f90_file, new_prob_str)
    end
end

function _build_path(base_path, ::Val{:dfmo_results}, problem_name, constraint_index)
	return joinpath(_build_subpath(base_path, problem_name, constraint_index), "dfmo_results.jld2")
end