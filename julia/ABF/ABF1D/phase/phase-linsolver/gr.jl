using Random
using BandedMatrices
using LinearAlgebra, SparseArrays
using Distributions
using Distributed
using Parameters
include("./phase.jl")

function compute_psi(; L::Integer = 1001,θ::F = 0.25, rng::AbstractRNG = Random.GLOBAL_RNG) where F <: AbstractFloat
    H = ham_sf_obc(L = L, θ = θ, rng = rng)
    H = Complex{F}.(BandedMatrix(H, (2, 2)));    
    Y = Array(@view H[:, 1])
    A = H[:, 2:end]
    X = A\Y
    pushfirst!(X, 1.)    
    return X
end

function eig_corr!(y, x, i, len_cor)
      @views y .+= abs.(x[i:(i+len_cor-1)] .- x[i])
end

function eig_corr(y, x, i, len_cor)
      y = zeros(Float64, len_cor)
      eig_corr!(y, x, i, len_cor)
end
"""
Correlation averaged over j and realization
"""
function cor_tmm(;θ::F = 0.25,  L::Integer = 10000, R::Integer = 10, rng::AbstractRNG = Random.GLOBAL_RNG) where F <: AbstractFloat
    cutoff_wf_1 = 1
    cutoff_wf_2 = L
    len_cutoff = cutoff_wf_2 - cutoff_wf_1 + 1
    len_r = len_cutoff÷2
    g_r = zeros(Float64, len_r)
    g_r_sq = zeros(Float64, len_r)
    
    for r in 1:R
        data = compute_psi(θ = θ, L = L, rng = rng)
        logpsi = Float64.(log.(abs.(data)))
        logpsi_sq = logpsi.^2
        for i in 1:len_r
            eig_corr!(g_r, logpsi, i, len_r)
            eig_corr!(g_r_sq, logpsi_sq, i, len_r)
        end
    end
    
     return g_r/R/len_r, g_r_sq/R/len_r
end

"""
(parallel version) Correlation averaged over j and realization. This is called whenever the input argument rng is an array of rngs.
"""
function cor_tmm(;θ::F = 0.25,  L::Integer = 10000, R::Integer = 10, rng::Vector{rngs} = [Random.GLOBAL_RNG]) where {F <: AbstractFloat, rngs <: AbstractRNG}
    cutoff_wf_1 = 1
    cutoff_wf_2 = L
    len_cutoff = cutoff_wf_2 - cutoff_wf_1 + 1
    len_r = len_cutoff÷2
    
    gr_sum = @distributed (+) for r in 1:R
        logpsi = Float64.(log.(abs.(compute_psi(θ = θ, L = L, rng = rng[myid()]))))
        g_r = zeros(Float64, len_r)
        for i in 1:len_r
            eig_corr!(g_r, logpsi, i, len_r)
        end
        g_r
    end
    
     return gr_sum/R/len_r 
end

@with_kw struct Params{F <: AbstractFloat}
    θ::F = 0.25
    R::Int64 = 2 
    L::Int64 = 10000
    seed::Int64 = 1234
end    



function read_config(dt::Dict)
    if haskey(dt, "vartype")
        s = dt["vartype"]
        F = eval(Meta.parse(s))
    else 
        F = Float64
    end
    return Params(θ = F(dt["th"]),  R = dt["R"], L = dt["L"], seed = dt["seed"])
end
