using Random
using LinearAlgebra
using DataFrames
include("abf2_ham.jl")
include("abf2_cores.jl")

@doc """
A struct where parameters are stored.
"""
struct Parameters{F<:AbstractFloat, N<:Int64}
    N::N # Number of unit cell.
    W::F
    V::F # If true, V = 10^V
    θ::F
    vl::F # Eigenvalues lower limits and upper limits. If no limit, vl = nothing, vu = nothing
    vu::F
    seed::N
    R::N
    core::N # 1: SF, on+phase, 2: SF, only phase, 3: FD, on+phase, 4: FD, only phase
end



@doc """
A struct where parameters are stored. This is the input parameters for the scanning.
    θ is in unit of π rad.
"""
struct Full_Parameters{F<:AbstractFloat,N<:Int64, B<:Bool}
    N::N # Number of unit cell.
    W::F
    V_min::F
    V_max::F
    V_num::N
    V_log::B # If true, V = 10^V
    θ_min::F
    θ_max::F
    θ_num::N
    vl::F # Eigenvalues lower limits and upper limits. If vl = vu, there is no limit
    vu::F
    seed::N
    R::N
    core::N # 1: SF, on+phase, 2: SF, only phase, 3: FD, on+phase, 4: FD, only phase
end

# Choose basis
function core_params(p::Parameters)
    if p.core == 1
        f = pn_sf_ph_on
    elseif p.core == 2
        f = pn_sf_ph
    elseif p.core == 3
        f = pn_fd_ph_on
    elseif p.core == 4
        f = pn_fd_ph
    end
    return f
end


## Full_Parameters -> Parameters
function full_to_one(p::Full_Parameters, i, j)
    V, θ = expand_params(p)
    Parameters(p.N, p.W, V[i], θ[j], p.vl, p.vu, p.seed, p.R, p.core)
end

@doc """
Read configuration JSON file, and construct a Parameters struct from it.
"""
function readconfig(config::Dict)
    p = Full_Parameters(
        config["N"],
        config["W"],
        config["delta_min"],
        config["delta_max"],
        config["delta_num"],
        config["delta_log"],
        config["theta_min"],
        config["theta_max"],
        config["theta_num"],
        config["vl"],
        config["vu"],
        config["seed"],
        config["R"],
        config["core"]
        )

    return p
end

@doc """
Expand domain parameters
"""
function expand_params(p::Full_Parameters)
    θ = collect(range(p.θ_min, p.θ_max, length = p.θ_num))
    V = collect(range(p.V_min, p.V_max, length = p.V_num))
    if p.V_log
        V = 10 .^V
    end
    return V, θ
end


@doc """
Scan obtain E and PN for R different realizations
"""
function scan_pn(p::Parameters)
    #-------------------Initialize output array------------------------------------#
    df_th = [DataFrame(E = Float64[], PN = Float64[]) for i in 1:Threads.nthreads()]
    #-------------------Prepare Hamiltonians and Unitaries-------------------------#
    H_FE, U = ham_FE(p.θ, p.N)
    #-------------------Scan setup------------------------------------------------------#
    rng = [MersenneTwister(p.seed+j) for j in 1:Threads.nthreads()]
    f = core_params(p)
    limit_ev = (p.vu == p.vl)
    #------------------- Scan -------------------#
    r_on = [Array{Float64}(undef, 2p.N) for j in 1:Threads.nthreads()]
    r_off = [Array{Float64}(undef, 12p.N) for j in 1:Threads.nthreads()]

    @Threads.threads for n = 1:p.R
        thid = Threads.threadid()
        r_on[thid] .= rand(rng[thid], 2p.N) .- 0.5
        r_off[thid] .= rand(rng[thid], 12p.N) .- 0.5

        D = spdiagm(0 => r_on[thid]) #diagonal disorder
        T = off_phase_dis2(H_FE, r_off[thid]) #off-diagonal phase disorder

        if limit_ev
            E, PN = f(H_FE, D, T, U, p.W, p.V)
        else
            E, PN = f(H_FE, D, T, U, p.W, p.V, p.vl, p.vu)
        end

        for i in 1:length(E)
            push!(df_th[thid], (pop!(E), pop!(PN)))
        end
    end
    df_out = reduce(vcat, df_th)

    return df_out
end
