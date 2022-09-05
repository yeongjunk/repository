module Typical 
export MFAParameters, prepare_MFA!, compute_fα, compute_α, compute_τ, compute_ταf
using Lattices
using Statistics
import LogExpFunctions: xlogx, xlogy
include("./box_counting.jl")

struct MFAParameters{Lat <: Lattice}
    ltc::Lat
    l::Vector{Int64}
    q::Vector{Float64}
    box_indices::Vector{Array{Int64, 2}}
    function MFAParameters(ltc::Lat, l, q) where {Lat <: Lattice}
        new{Lat}(ltc, l, q, Vector{Array{Int64, 2}}[]) 
    end
end

function prepare_MFA!(params::MFAParameters)
    for li in params.l
       idx_li = box_indices(params.ltc, li)
       push!(params.box_indices, idx_li)
    end
end

function compute_fα(p, L, l, q)
    μ = p.^q
    Z_q = sum(μ)
    for i in 1:size(p, 2)
        μ[:, i] ./= Z_q[i]
    end
    return (1/log(l/L))*sum(xlogx.(μ))
end

function compute_fα(params::MFAParameters, eigvect::Vector)
    p = abs2.(eigvect)
    p_coarse = [box_coarse(p, params.box_indices[i]) for i in 1:length(params.l)]
    fα = [compute_fα(p_coarse[i], L, params.l[i], params.q[j]) for i in 1:length(params.l), j in 1:length(params.q)]
    return fα
end

function compute_α(p, L, l, q)
    μ = p.^q
    Z_q = sum(μ)
    for i in 1:size(p, 2)
        μ[:, i] ./= Z_q[i]
    end
    return (1/log(l/L))*sum(xlogy.(μ, p))
end

function compute_α(params::MFAParameters, eigvect::Vector)
    p = abs2.(eigvect)
    p_coarse = [box_coarse(p, params.box_indices[i]) for i in 1:length(params.l)]
    α = [compute_α(p_coarse[i], L, params.l[i], params.q[j]) for i in 1:length(params.l), j in 1:length(params.q)]
    return α
end

function compute_τ(p::Vector, L, l, q)
    gipr = sum(x -> x^q, p)
    return log(gipr)/log(l/L) 
end

function compute_τ(params::MFAParameters, eigvect::Vector)
    p = abs2.(eigvect)
    p_coarse = [box_coarse(p, params.box_indices[i]) for i in 1:length(params.l)]
    τ = [compute_τ(p_coarse[i], L, params.l[i], params.q[j]) for i in 1:length(params.l), j in 1:length(params.q)]

    return τ
end

function compute_ταf(p::Vector, L, l, q; test=true) where F
    μ = p.^q
    Z_q = sum(μ)
    for i in 1:size(p, 2)
        μ[:, i] ./= Z_q[i]
    end
    α = sum(xlogy.(μ, p))/log(l/L)
    f_α = (1/log(l/L))*sum(xlogx.(μ))/log(l/L)

    if test 
        τ = α*q - f_α 
    else
        τ =log(Z_q)/log(l/L)
    end

    return log(Z_q)/log(l/L), (1/log(l/L))*sum(xlogy.(μ, p)), (1/log(l/L))*sum(xlogx.(μ))

end

function compute_ταf(params::MFAParameters, eigvect::Array{F, 1}) where F
    p = abs2.(eigvect)
    p_coarse = [box_coarse(p, params.box_indices[i]) for i in 1:length(params.l)]
    full = [compute_ταf(p_coarse[i], params.ltc.N, params.l[i], params.q[j]) for j in 1:length(params.q), i in 1:length(params.l)]
    τ = getindex.(full, 1)
    α = getindex.(full, 2)
    fα = getindex.(full, 3)
    return τ, α, fα 
end

end # module
