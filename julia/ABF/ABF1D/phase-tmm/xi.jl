using Random
using LinearAlgebra
using ProgressBars
using LaTeXStrings
using CSV, DataFrames
using Distributions
using Parameters
using StatsBase
@doc """
t_1: nearest neighbor, t_2: next nearest neighbor hopping for scale free model of RMF disorder.
p is 5 by 2 array containing 10 random variables.
"""
function hoppings(p, s, c, s4,c4, s2c2)
    t_1 = 2im*(p[1, 1]*s4 - p[1, 1]*c4 + p[2, 1]*s2c2 + p[2, 2]*s2c2 
        + p[3, 1]*c4 - p[3, 2]*s2c2 - p[4, 1]*s2c2 + p[4, 2]*s4
        - p[5, 1]*c4 - p[5, 2]*s4)*s2c2
    t_2 = 2im*(-p[2, 1] + p[3, 1] + p[4, 1] - p[5, 1])*s4*c4
    return t_1, t_2
end

@doc """
Overwrite transfer matrix. 
"""
function tm!(TM, E, t_1, t_2, t_m1, t_m2, tt_1, tt_2, tt_m1, tt_m2)
    TM[1, :] .= [(E/tt_2); (-tt_1/tt_2); (-tt_m2/tt_2); (-tt_m1/tt_2)]
    TM[2, :] .= [(-E*t_1/tt_2-t_m1)/t_2; (E + t_1*tt_1/tt_2)/t_2; (t_1*tt_m2/tt_2)/t_2; (t_1*tt_m1/tt_2 - t_m2)/t_2]
    TM[3:4, 1:2] .= I(2)
    TM[3:4, 3:4] .= 0
end

@doc """
Generate transfer matrix.
"""
function tm(E, t_1, t_2, t_m1, t_m2, tt_1, tt_2, tt_m1, tt_m2, vartype)
    TM = Array{Complex{vartype}}(undef, 4, 4)
    tm!(TM, E, t_1, t_2, t_m1, t_m2, tt_1, tt_2, tt_m1, tt_m2)

    return TM
end

@doc """
for T_j, the transfer vector is given by (ψ_j-1, ψ_j, ψ_j-3, ψ_j-2). j should increase by 2 in each step
"""
function compute_xi(;θ::Real = 0.25, E::Real = 0.1 , N::Integer = 100000, rng::AbstractRNG = Random.GLOBAL_RNG)
    vartype = typeof(θ)
    dist = Uniform(vartype(-0.5), vartype(0.5))
    # pre-calculation of sines and cosines for performance optimization
    s, c = sincospi(θ)
    s4 = s^4; c4 = c^4; s2c2 = s^2*c^2
    # Initial hoppings
    p = rand(rng, dist, (5, 2))
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
    ψⱼ = Complex{vartype}.([ψ₂ ψ₂_2; ψ₃ ψ₃_2; ψ₁ ψ₁_2; ψ₀ ψ₀_2])
    ΣlnR = vartype(0.);
    T = Array{Complex{vartype}}(undef, 4, 4)
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
        @views ψⱼ[1:2,2] .-= dot(ψⱼ[1:2,1], ψⱼ[1:2,2])/dot(ψⱼ[1:2,1], ψⱼ[1:2,1])*ψⱼ[1:2,1] # Gram-Schmidt orthogonalization
        ΣlnR += log(abs(ψⱼ[1, 2]/ψⱼ[3, 2]))
        normalize!(@view ψⱼ[:, 2])
        normalize!(@view ψⱼ[:, 1])
    end
    ξ = N / ΣlnR
    return ξ
end

@with_kw struct Params
    θ::Float64 = 0.25
    E::Vector{Float64} = [0.1]
    R::Int64 = 1
    N::Vector{Int64} = [10000]
    seed::Int64 = 1234
end    

function read_config(dt::Dict)
    return Params(dt["th"], dt["E"], dt["R"], dt["N"], dt["seed"])
end

function scan_xi(p::Params)
    @unpack θ, E, R, N, seed = p
    rng = [MersenneTwister(seed + i) for i in 1:Threads.nthreads()]
    xi = Array{Float64}(undef, length(E), R) 
    for i in tqdm(1:length(E))
        @Threads.threads for r in 1:R
            xi[i, r] = compute_xi(θ = θ, E = E[i], N = N[i], rng = rng[r])
        end
    end
    return mean(xi, dims = 2)
end
