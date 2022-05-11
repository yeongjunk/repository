using Random
using LinearAlgebra, SparseArrays
using ProgressBars
using CSV, DataFrames
using Distributions
using Parameters
using StatsBase
@doc """
t_1: nearest neighbor, t_2: next nearest neighbor hopping for scale free model of RMF disorder.
p is 5 by 2 array containing 10 random variables.
"""
function hoppings(p::AbstractArray{T}, s::T, c::T, s4::T,c4::T, s2c2::T) where T <: Real
    t_1 = 2im*(p[1, 1]*s4 - p[1, 1]*c4 + p[2, 1]*s2c2 + p[2, 2]*s2c2 
        + p[3, 1]*c4 - p[3, 2]*s2c2 - p[4, 1]*s2c2 + p[4, 2]*s4
        - p[5, 1]*c4 - p[5, 2]*s4)*s2c2
    t_2 = 2im*(p[2, 1] - p[3, 1] - p[4, 1] + p[5, 1])*s4*c4
    return t_1, t_2
end

function tm!(TM::AbstractArray{T}, E::F, t_1::T, t_2::T, t_m1::T, t_m2::T, tt_1::T, tt_2::T, tt_m1::T, tt_m2::T) where {T <: Number, F <: Number} 
    TM[1, 1] = E/tt_2
    TM[1, 2] = -tt_1/tt_2
    TM[1, 3] = -tt_m2/tt_2
    TM[1, 4] = -tt_m1/tt_2
    TM[2, 1] = (-E*t_1/tt_2-t_m1)/t_2
    TM[2, 2] = (E + t_1*tt_1/tt_2)/t_2
    TM[2, 3] = (t_1*tt_m2/tt_2)/t_2
    TM[2, 4] = (t_1*tt_m1/tt_2 - t_m2)/t_2
    TM[3, 1] = one(T) 
    TM[4, 2] = one(T)
end

@doc """
for T_j, the transfer vector is given by (ψ_j-1, ψ_j, ψ_j-3, ψ_j-2). j should increase by 2 in each step
"""
function compute_psi(;θ::F = 0.25, E::F = 0.1 , N::Integer = 100000, q = 2, rng::AbstractRNG = Random.GLOBAL_RNG) where F <: Real
    dist = Uniform(F(-0.5), F(0.5)) # Uniform distribution [-0.5, 0.5]
    s = sinpi(θ)
    c = cospi(θ)
    s4 = s^4; c4 = c^4; s2c2 = s^2*c^2 # pre-calculation of sines and cosines for performance optimization
    # Initial hoppings
    p = Array{F}(undef, 5, 2)
    rand!(rng, dist, p)
    t₁ₙ₋₃, t₂ₙ₋₃ = hoppings(p, s, c, s4, c4, s2c2)
    p[:, 2] = @view p[:, 1]
    rand!(rng, dist, @view p[:, 1])
    t₁ₙ₋₂, t₂ₙ₋₂ = hoppings(p, s, c, s4, c4, s2c2)
    p[:, 2] = @view p[:, 1]
    rand!(rng, dist, @view p[:, 1])
    t₁ₙ₋₁, t₂ₙ₋₁ = hoppings(p, s, c, s4, c4, s2c2)
    p[:, 2] .= @view p[:, 1]
    rand!(rng, dist, @view p[:, 1])
    t₁ₙ, t₂ₙ = hoppings(p, s, c, s4, c4, s2c2)
    p[:, 2] .= @view p[:, 1]
    rand!(rng, dist, @view p[:, 1])
    
    # Initialization of vector. ψ₀,  ψ₁ cannot be known. and set to random. 
    ψ₀ = 0. 
    ψ₁ = 1. 
    ψ₂ = (t₂ₙ₋₃*ψ₀ + t₁ₙ₋₂*ψ₁)/t₂ₙ₋₁  
    ψ₃ = (t₂ₙ₋₂*ψ₁ + t₁ₙ₋₁*ψ₂)/t₂ₙ  
    ψ₀_2 = 1. 
    ψ₁_2 = 0. 
    ψ₂_2 = (t₂ₙ₋₃*ψ₀_2 + t₁ₙ₋₂*ψ₁_2)/t₂ₙ₋₁  
    ψ₃_2 = (t₂ₙ₋₂*ψ₁_2 + t₁ₙ₋₁*ψ₂_2)/t₂ₙ  
    ψⱼ = Complex{F}.([ψ₂ ψ₂_2; ψ₃ ψ₃_2; ψ₁ ψ₁_2; ψ₀ ψ₀_2]) 
   
    
    # Initialization of vector. (ψ₀,  ψ₁) cannot be known, set to orthogonal vectors (0, 1) and (1, 0)
    T = zeros(Complex{F}, 4, 4)
    data = Complex{F}[]


    
    for j in 1:N
        tm!(T, E, t₁ₙ, t₂ₙ, conj(t₁ₙ₋₁), conj(t₂ₙ₋₂), t₁ₙ₋₁, t₂ₙ₋₁, conj(t₁ₙ₋₂), conj(t₂ₙ₋₃))
        ψⱼ .= T*ψⱼ
        t₁ₙ₋₃, t₂ₙ₋₃, t₁ₙ₋₂, t₂ₙ₋₂ = t₁ₙ₋₁, t₂ₙ₋₁, t₁ₙ, t₂ₙ
        p[:, 2] .= @view p[:, 1]
        rand!(rng, dist, @view p[:, 1])
        t₁ₙ₋₁, t₂ₙ₋₁ = hoppings(p, s, c, s4, c4, s2c2)
        p[:, 2] .= @view p[:, 1]
        rand!(rng, dist, @view p[:, 1])
        t₁ₙ, t₂ₙ = hoppings(p, s, c, s4, c4, s2c2)
        if j%q == 0
            ψ_new, R = qr!(ψⱼ)
            ψⱼ .= Matrix{Complex{F}}(ψ_new)
            @views ψⱼ[:, 1] .*= R[1, 1]
            @views ψⱼ[:, 2] .*= R[2, 2]

            push!(data, ψⱼ[1, 2])
            push!(data, ψⱼ[2, 2])
        end
    end

    return data
end

function eig_corr!(y, x, i, len_cor)
      @views y .+= abs.(x[i:(i+len_cor-1)] .- x[i])
end

"""
Correlation averaged over j and realization
"""
function cor_tmm(;E = 0., θ = 0.25, δ = 0., N = 10000, q = 4, R = 10, rng = Random.GLOBAL_RNG)
    cutoff_wf_1 = 1
    cutoff_wf_2 = N
    len_cutoff = cutoff_wf_2 - cutoff_wf_1 + 1
    len_r = len_cutoff÷2
    g_r = zeros(Float64, len_r)
    g_r_sq = zeros(Float64, len_r)
    
    for r in 1:R
        data = compute_psi(θ = θ, E = E, N = N, q = q, rng = rng)
        logpsi = Float64.(log.(abs.(data)))
        logpsi_sq = logpsi.^2
        for i in 1:len_r
            eig_corr!(g_r, logpsi, i, len_r)
            eig_corr!(g_r_sq, logpsi_sq, i, len_r)
        end
    end
    
     return g_r/R/len_r, g_r_sq/R/len_r
end

@with_kw struct Params{F <: AbstractFloat}
    θ::F = 0.25
    E::Vector{F} = [0.1]
    R::Int64 = 2 
    N::Int64 = 10000
    q::Int64 = 2
    seed::Int64 = 1234
end    



function read_config(dt::Dict)
    if haskey(dt, "vartype")
        s = dt["vartype"]
        F = eval(Meta.parse(s))
    else 
        F = Float64
    end
    if !haskey(dt, "q")
       dt["q"] = 2
    end 
    return Params(θ = F(dt["th"]), E = Vector{F}(dt["E"]), R = dt["R"], N = dt["N"], seed = dt["seed"], q = dt["q"])
end
