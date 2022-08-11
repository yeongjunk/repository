using Random
using BandedMatrices
using LinearAlgebra, SparseArrays
using Distributions
using Distributed
using Parameters
include("./phase.jl")

struct LatticeWave{F}
    lattice::Lattice2D
    ψ::Vector{F} 
    function LatticeWave(ltc::Lattice2D, psi::Vector{F}) where F
        new{F}(ltc, psi) 
    end
end
function compute_psi(; L::Integer = 1001,θ::F = 0.25, rng::AbstractRNG = Random.GLOBAL_RNG) where F <: AbstractFloat
    H = ham_sf_obc(L = L, θ = θ, rng = rng)
    H = Complex{F}.(BandedMatrix(H, (2L+2, 2L+2)));
    Y = Array(@view H[:, 1])
    A = H[:, 2:end]
    X = A\Y
    pushfirst!(X, 1.)    
    return X
end

"""
Take the log|ψ| as input. compute the g(r) at an index (m, n)
"""
function eig_corr(lw::LatticeWave, (m,n); direction = :x)
    @assert lw.lattice.U == 1
    ltc = lw.lattice
    i = index(ltc, (m, n, 1))
    if direction == :x
        f = [index(ltc, (m+j, n, 1)) for j in 1:ltc.M-m]
    elseif direction == :y
        f = [index(ltc, (m, n+j, 1)) for j in 1:ltc.N-n]
    else
        error("direction has to be either :x and :y")
    end

    return vec(abs.(lw.ψ[f] .- lw.ψ[i]))
end
"""
Take log|ψ| as input. Compute the g(r) averaged over all points (m, n)
"""
function eig_corr_full(lw::LatticeWave{F}) where F
    @assert lw.lattice.U == 1
    @assert lw.lattice.M == lw.lattice.N
    ltc = lw.lattice
    grs = Vector{Float64}[]
    for i in 1:length(lw.ψ)
        (m,n,_) = site(ltc, i)
        gr_mn_x = eig_corr(lw, (m,n), direction = :x) 
        gr_mn_y = eig_corr(lw, (m,n), direction = :x) 
      
        gr_mn_x = [gr_mn_x; fill(F(NaN), ltc.M-length(gr_mn_x))]    
        gr_mn_y = [gr_mn_y; fill(F(NaN), ltc.M-length(gr_mn_y))]    
        push!(grs, gr_mn_x)
        push!(grs, gr_mn_y)
    end
    grs = reduce(hcat, grs) 
    grs_mean = Array{Float64}(undef, size(grs, 1))
    for i in 1:size(grs, 1)
        @views grs_mean[i] = mean(filter(!isnan, grs[i, :]))
    end
    return grs_mean
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
