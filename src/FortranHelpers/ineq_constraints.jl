
const constraint0_f90 = """
subroutine fconstriq(n,m,x,ciq)
	implicit none
	integer, intent(in) :: n,m
	real*8, intent(in), dimension(n) :: x
	real*8, intent(out), dimension(m) :: ciq

	return

end subroutine fconstriq
"""
fortran_constraint_string(::Val{0})=constraint0_f90
julia_constraint_func(::Val{0})=(ciq,x) -> nothing
num_vars_to_num_constraints(::Val, n_vars) = 0

"""
For ``x ∈ ℝ^n``, compute in-place
```math
	g_j = (3-2x_{j+1})x_{j+1} - x_j - 2x_{j+2} + 1,
```
for all ``j=1, …, m`` with ``m=n-2``.
"""
@views function constraint1!(ciq, x)
	@. ciq = (3 - 2*x[2:end-1])*x[2:end-1] - x[1:end-2] - 2*x[3:end] + 1
	return nothing
end
const constraint1_f90 = """
subroutine fconstriq(n, m, x, ciq)
	implicit none
    integer, intent(in) :: n, m
    real*8, intent(in), dimension(n) :: x
	real*8, intent(out), dimension(m) :: ciq
    
    integer :: j

    do j = 1, m
        ciq(j) = (3 - 2*x(j+1))*x(j+1) - x(j) - 2*x(j+2) + 1.0D+00
    end do

end subroutine fconstriq
"""
julia_constraint_func(::Val{1})=constraint1!
num_vars_to_num_constraints(::Val{1}, n_vars) = n_vars - 2
fortran_constraint_string(::Val{1})=constraint1_f90

"""
For ``x ∈ ℝ^n``, compute in-place
```math
	g_j = (3-2x_{j+1})x_{j+1} - x_j - 2x_{j+2} + 2.5,
```
for all ``j=1, …, m`` with ``m=n-2``.
"""
@views function constraint2!(ciq, x)
	@. ciq = (3 - 2*x[2:end-1])*x[2:end-1] - x[1:end-2] - 2*x[3:end] + 2.5
	return nothing
end
const constraint2_f90 = """
subroutine fconstriq(n, m, x, ciq)
	implicit none
    integer, intent(in) :: n, m
    real*8, intent(in), dimension(n) :: x
	real*8, intent(out), dimension(m) :: ciq

    integer :: j

    do j = 1, m
        ciq(j) = (3 - 2*x(j+1))*x(j+1) - x(j) - 2*x(j+2) + 2.5D+00
    end do

end subroutine fconstriq
"""
julia_constraint_func(::Val{2})=constraint2!
num_vars_to_num_constraints(::Val{2}, n_vars) = n_vars - 2
fortran_constraint_string(::Val{2})=constraint2_f90

"""
For ``x ∈ ℝ^n``, compute in-place
```math
	g_j = x_j^2 + x_{j+1}^2 + x_jx_{j+1} - 2x_j - 2 x_{j+1} + 1
```
for all ``j=1, …, m`` with ``m=n-1``.
"""
@views function constraint3!(ciq, x)
	@. ciq = x[1:end-1]*(x[1:end-1] + x[2:end] - 2) + x[2:end]*(x[2:end] - 2) + 1
 #   @. ciq = x[1:end-1]^2 + x[2:end]^2 + x[1:end-1]*x[2:end] - 2*x[1:end-1] - 2*x[2:end] + 1
	return nothing
end
const constraint3_f90 = """
subroutine fconstriq(n, m, x, ciq)
	implicit none
    integer, intent(in) :: n, m
    real*8, intent(in), dimension(n) :: x
	real*8, intent(out), dimension(m) :: ciq

    integer :: j

    do j = 1, m
        ciq(j) = x(j)*(x(j) + x(j+1) - 2) + x(j+1)*(x(j+1) - 2) + 1.0D+00
    end do

end subroutine fconstriq
"""
julia_constraint_func(::Val{3})=constraint3!
num_vars_to_num_constraints(::Val{3}, n_vars) = n_vars - 1
fortran_constraint_string(::Val{3})=constraint3_f90


"""
For ``x ∈ ℝ^n``, compute in-place
```math
	g_j = x_j^2 + x_{j+1}^2 + x_jx_{j+1} - 1
```
for all ``j=1, …, m`` with ``m=n-1``.
"""
@views function constraint4!(ciq, x)
	@. ciq = x[1:end-1]^2 + x[2:end]^2 + x[1:end-1]*x[2:end] - 1
	return nothing
end
const constraint4_f90 = """
subroutine fconstriq(n, m, x, ciq)
	implicit none
    integer, intent(in) :: n, m
    real*8, intent(in), dimension(n) :: x
	real*8, intent(out), dimension(m) :: ciq

    integer :: j

    do j = 1, m
        ciq(j) = x(j)**2 + x(j+1)**2 + x(j)*x(j+1) - 1.0D+00
    end do

end subroutine fconstriq
"""
julia_constraint_func(::Val{4})=constraint4!
num_vars_to_num_constraints(::Val{4}, n_vars) = n_vars - 1
fortran_constraint_string(::Val{4})=constraint4_f90

"""
For ``x ∈ ℝ^n``, compute in-place
```math
	g_j = (3 - 0.5x_{j+1})x_{j+1} - x_j - 2x_{j+2} + 1
```
for all ``j=1, …, m`` with ``m=n-2``.
"""
@views function constraint5!(ciq, x)
	@. ciq = (3 - 0.5*x[2:end-1])*x[2:end-1] - x[1:end-2] - 2*x[3:end] + 1
	return nothing
end
const constraint5_f90 = """
subroutine fconstriq(n, m, x, ciq)
	implicit none
    integer, intent(in) :: n, m
    real*8, intent(in), dimension(n) :: x
	real*8, intent(out), dimension(m) :: ciq

    integer :: j

    do j = 1, m
        ciq(j) = (3 - 0.5*x(j+1))*x(j+1) - x(j) - 2*x(j+2) + 1.0D+00
    end do

end subroutine fconstriq
"""
julia_constraint_func(::Val{5})=constraint5!
num_vars_to_num_constraints(::Val{5}, n_vars) = n_vars - 2
fortran_constraint_string(::Val{5})=constraint5_f90

"""
For ``x ∈ ℝ^n``, compute in-place
```math
	g_1 = \\sum_{j=1}^{n-2}((3 - 0.5x_{j+1})x_{j+1} - x_j - 2x_{j+2} + 1)
```
"""
@views function constraint6!(ciq, x)
	ciq[end] = 0
	for j=1:length(x)-2
		ciq[end] += (3 - 0.5*x[j+1])*x[j+1] - x[j] - 2*x[j+2] + 1
	end
	return nothing
end
const constraint6_f90 = """
subroutine fconstriq(n, m, x, ciq)
	implicit none
    integer, intent(in) :: n, m
    real*8, intent(in), dimension(n) :: x
	real*8, intent(out), dimension(m) :: ciq

    integer :: j

    ciq(1) = 0.0D+00
    do j = 1, n-2
        ciq(1) = ciq(1) + (3 - 0.5*x(j+1))*x(j+1) - x(j) - 2*x(j+2) + 1.0D+00
    end do

end subroutine fconstriq
"""
julia_constraint_func(::Val{6})=constraint6!
num_vars_to_num_constraints(::Val{6}, n_vars) = 1
fortran_constraint_string(::Val{6})=constraint6_f90

const f90_strings = (
    constraint0_f90,
    constraint1_f90,
    constraint2_f90,
    constraint3_f90,
    constraint4_f90,
    constraint5_f90,
    constraint6_f90
)
