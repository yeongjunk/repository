using Distributions
using LinearAlgebra
using Random
using Statistics
using ROAG
using StatsBase
using Parameters

"""
Symmetric GOE matrix.
"""
function GOE(;rng = Random.GLOBAL_RNG, σ = 1., q = 10)
    M = rand(rng, Normal(0, σ), q, q)
    return Hermitian((M + M') / 2)
end

wigner_dyson(s) = pi/2*s*exp(-pi*s^2/4)
poisson(s) = exp(-s)

"""
Floquet Rosenzweig-Porter uniform
"""
function FRP_uniform(;rng = Random.GLOBAL_RNG, q = 2, σ = 1., γ = 1., g = 1)
    Λ = Diagonal(rand(rng, Uniform(-σ, σ), q))
    G = GOE(σ = 1., rng = rng, q = q)
    return exp(-im*Λ)*exp(-im*(g/(q^γ)*G))
end

"""
Floquet Rosenzweig-Porter Gaussian
"""
function FRP_normal(;rng = Random.GLOBAL_RNG, q = 2, σ = 1., γ = 1., g = 1)
    Λ = Diagonal(rand(rng, Normal(0, σ), q))
    G = GOE(σ = 1., rng = rng, q = q)
    return exp(-im*Λ)*exp(-im*(g/(q^γ)*G))
end

"""
Obtain mean of ROAG over R different realizations of q X q FRP matrix(Λ uniform distribution)
"""
function roag_FRP_normal(;rng = Random.GLOBAL_RNG, q = 5, σ = 1., γ = 1., g = 1., R = 1000)
    roag_mean = Float64[]
    for i in 1:R
        exp_iE = eigvals!(FRP_normal(rng = rng, q = q, σ = σ, γ = γ, g = g))
        E = sort(angle.(exp_iE))
        roag!(E)
        push!(roag_mean, mean(E))
    end
    return mean(roag_mean), std(roag_mean)/sqrt(R)
end


function roag_FRP_uniform(;rng = Random.GLOBAL_RNG, q = 5, σ = 1., γ = 1., g = 1., R = 1000)
    roag_mean = Float64[]
    for i in 1:R
        exp_iE = eigvals!(FRP_uniform(rng = rng, q = q, σ = σ, γ = γ, g = g))
        E = sort(angle.(exp_iE))
        roag!(E)
        push!(roag_mean, mean(E))
    end
    return mean(roag_mean), std(roag_mean)/sqrt(R)
end

##---------- Scan functions----------##

@with_kw struct ScanParameters{A1 <: AbstractArray, A2 <: AbstractArray, A3 <: AbstractArray, F}
   γ::A1 = [0.0, 0.1, 0.3, 0.4, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
   q::A2 = [10, 15, 20]
   R::A3 = [10000, 10000, 5000]
   σ::F = 1.0
   g::F = 1.0
   seed::Int64 = 1234
   dist_diag::Bool = true # true: normal, false: uniform
end

function dict_to_params(dict)
    return ScanParameters(γ = dict["gamma"], q = dict["q"], R = dict["R"], σ = dict["sigma"], g = dict["g"], seed = dict["seed"], dist_diag = dict["dist_diag"])
end

function params_to_dict(p::ScanParameters)
    @unpack γ, q, R, σ, g, seed = p 
    return Dict("gamma" => γ, "q" => q, "R" => R, "sigma" => σ, "g" => g, "seed" => seed)
end

function roag_FRP_scan(p::ScanParameters)
    @unpack γ, q, R, σ, g, seed, dist_diag = p
    rng = MersenneTwister(seed)
    R_mean = Array{Float64}(undef, length(γ), length(q))
    R_ste = similar(R_mean)
    for i in 1:length(q)
        println("start q = ", q[i])
        @time for j in 1:length(γ)
            if dist_diag == true
                R_mean[j ,i], R_ste[j, i] = roag_FRP_normal(rng = rng, q = q[i], σ = σ, g = g, γ = γ[j], R = R[i])
            else
                R_mean[j ,i], R_ste[j, i] = roag_FRP_uniform(rng = rng, q = q[i], σ = σ, g = g, γ = γ[j], R = R[i])
            end
        end
    end

    dct = params_to_dict(p)
    dct["R_mean"] = R_mean
    dct["R_ste"] = R_ste

    return dct
end

