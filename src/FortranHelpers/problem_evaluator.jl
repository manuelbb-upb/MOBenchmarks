
function ccall_dims(ptr::Ptr{Cvoid})
    num_vars = Ref(0) :: Ref{Int}				# `n` in Fortran
    num_nl_ineq_constrs = Ref(0) :: Ref{Int} 	# `m` in Fortran
    num_objfs = Ref(0) :: Ref{Int}				# `q` in Fortran

    ccall(
        ptr, 
        Nothing,
        (Ref{Int}, Ref{Int}, Ref{Int}),
        num_vars, num_nl_ineq_constrs, num_objfs
    )

    return num_vars[], num_nl_ineq_constrs[], num_objfs[]
end

function ccall_startp(ptr::Ptr{Cvoid}, n_vars)
	x0 = zeros(Float64, n_vars)
	ccall(
        ptr,
		Nothing,
		(Ref{Int}, Ptr{Float64}),
		n_vars, x0	# does not seem to matter that we pass an Int instead of Ref{Int}
	)
	return x0
end#function get_startp

function ccall_bounds(ptr::Ptr{Cvoid}, n_vars)
	lb = zeros(Float64, n_vars)
	ub = copy(lb)
	ccall(
		ptr,
		Nothing,
		(Ref{Int}, Ptr{Float64}, Ptr{Float64}),
	    n_vars, lb, ub
	)
	return lb, ub
end

struct ProblemEvaluator
    dl :: Ptr{Cvoid}
    setdim_ptr :: Union{Ptr{Cvoid}, Nothing}
    startp_ptr :: Union{Ptr{Cvoid}, Nothing}
    setbounds_ptr :: Union{Ptr{Cvoid}, Nothing}
    functs_ptr :: Union{Ptr{Cvoid}, Nothing}
    fconstriq_ptr :: Union{Ptr{Cvoid}, Nothing}

    n_vars :: Int
    n_constrs :: Int 
    n_objfs :: Int
    
    counter_objfs :: Base.RefValue{Int}
    counter_constrs :: Base.RefValue{Int}

    x0 :: Vector{Float64}
    
    lb :: Vector{Float64}
    ub :: Vector{Float64}

    c_ineq :: Vector{Float64}
end

function ProblemEvaluator(sd_path)
    dl = Libdl.dlopen(sd_path)
    setdim_ptr = Libdl.dlsym(dl, :setdim; throw_error=false)
    startp_ptr = Libdl.dlsym(dl, :startp; throw_error=false)
    setbounds_ptr = Libdl.dlsym(dl, :setbounds; throw_error=false)
    functs_ptr = Libdl.dlsym(dl, :functs; throw_error=false)
    fconstriq_ptr = Libdl.dlsym(dl, :fconstriq; throw_error=false)
    n_vars, n_constrs, n_objfs = if isnothing(setdim_ptr)
        (0, 0, 0)
    else
        ccall_dims(setdim_ptr)
    end
    x0 = isnothing(startp_ptr) ? zeros(n_vars) : ccall_startp(startp_ptr, n_vars)

    lb, ub = if isnothing(setbounds_ptr) 
        (fill(-Inf, n_vars), fill(Inf, n_vars))
    else
        ccall_bounds(setbounds_ptr, n_vars)
    end

    c_ineq = zeros(n_constrs)
    return ProblemEvaluator(
        dl, setdim_ptr, startp_ptr, setbounds_ptr, functs_ptr, fconstriq_ptr,
        n_vars, n_constrs, n_objfs, Ref(0), Ref(0), x0, lb, ub, c_ineq
    )
end

function fortran_wrap(func, sd_path)
    ret_val = nothing
    p = ProblemEvaluator(sd_path)
    try 
        ret_val = func(p)
    finally
        Libdl.dlclose(p.dl)
    end
    p = nothing
    GC.safepoint()
    #GC.gc(false)
    return ret_val
end

function dims(ev::ProblemEvaluator)
    if !isnothing(ev.setdim_ptr)
        return ccall_dims(ev.setdim_ptr)
    end
    return nothing
end
function startp(ev::ProblemEvaluator)
    if !isnothing(ev.startp_ptr)
        return ccall_startp(ev.startp_ptr, ev.n_vars)
    end 
    return nothing
end

bounds(ev::ProblemEvaluator)=ccall_bounds(ev.setbounds_ptr, ev.n_vars)
function inplace_objectives(ev::ProblemEvaluator)
    n_vars = ev.n_vars
    n_objfs = ev.n_objfs
    ev.counter_objfs[] = 0
    return function(y, x)
        ev.counter_objfs[] += 1
        ccall(
            ev.functs_ptr,
            Nothing,
            (Ref{Int}, Ptr{Float64}, Ref{Int}, Ptr{Float64}),
            n_vars, x, n_objfs, y
        )
        return nothing
    end
end
function inplace_constraints(ev::ProblemEvaluator)
    n_vars = ev.n_vars
    n_constrs = ev.n_constrs
    ev.counter_constrs[] = 0
    return function(y, x)
        ev.counter_constrs[] += 1
        ccall(
            ev.fconstriq_ptr,
            Nothing,
            (Ref{Int}, Ref{Int}, Ptr{Float64}, Ptr{Float64}),
            n_vars, n_constrs, x, y
        )
        return nothing
    end
end

function l1_penalized_objectives(ev::ProblemEvaluator; eps=0.1)
    objfs! = inplace_objectives(ev)
    constrs! = inplace_constraints(ev)
    n_constrs = ev.n_constrs
    cvec = zeros(n_constrs)
    function (y, x)
        objfs!(y, x)
        constrs!(cvec, x)
        viol = _eps_term(cvec, eps)
        y .+= viol
        return nothing
    end
end

function _eps_term(cx, eps)
    pen = 0
    for (i,ci) in enumerate(cx)
        ci <= 0 && continue
        pen += (ci / _penalty_eps(eps, i))
    end
    return pen
end

_penalty_eps(eps::Number, i)=eps
_penalty_eps(eps::AbstractVector, i) = eps[i]

# f(x) + eps * sum( max.(0, c(x)).^(1/gamma) )
function lu_penalized_objectives(ev::ProblemEvaluator; eps=0.01, gamma=2)
    @assert gamma > 1
    objfs! = inplace_objectives(ev)
    constrs! = inplace_constraints(ev)
    n_constrs = ev.n_constrs
    cvec = zeros(n_constrs)
    function (y, x)
        objfs!(y, x)
        constrs!(cvec, x)
        viol = sum(eps .* (max.(0, cvec).^(1/gamma)))
        y .+= viol
        return nothing
    end
end

