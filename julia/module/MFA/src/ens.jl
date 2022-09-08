module Ensemble

export MFAParameters, prepare_MFA!, compute_gipr, compute_gipr_2, compute_τ, compute_ταf

using Lattices
using Statistics
import LogExpFunctions: xlogx, xlogy
include("./box_counting.jl")

dropmean(A; dims=:) = dropdims(mean(A; dims=dims); dims=dims)

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

function compute_gipr(params::MFAParameters, eigvect::AbstractArray{F, 1}) where F
    p = abs2.(eigvect)
    p_coarse = [box_coarse(p, params.box_indices[i]) for i in 1:length(params.l)]
    gipr = Array{Float64}(undef,  length(params.q), length(params.l))
    for i in 1:length(params.q), j in 1:length(params.l)
        gipr[i, j] = sum(x -> x^params.q[i], p_coarse[j])
    end
    return gipr
end

function compute_gipr_2(params::MFAParameters, eigvect::AbstractArray{F, 1}) where F
    @views p_coarse = [box_coarse(eigvect, params.box_indices[i], density = false) for i in 1:length(params.l)]

    gipr = zeros(Float64, length(params.q), length(params.l))
    μqlnμ= zeros(Float64, length(params.q), length(params.l))
    for j in 1:length(params.l), i in 1:length(params.q)
        for k in 1:length(p_coarse[j])
            p_q = p_coarse[j][k]^params.q[i] 
            gipr[i, j]  += p_q 
            μqlnμ[i, j] += xlogy(p_q, p_coarse[j][k])
        end
    end
    return gipr, μqlnμ
end

function compute_τ(params::MFAParameters, gipr::Array{F, 3}) where F
    gipr_mean = dropmean(gipr, dims=1)
    τ = similar(gipr_mean)
    for i in 1:size(τ, 2)
        τ[:, i] = log.(@view gipr_mean[:, i]) ./ log.(params.l[i]/params.ltc.N)
    end
    return τ
end

function compute_τ(params::MFAParameters, eigvects::AbstractArray{F, 2}) where F
    gipr = Array{Float64}(undef, size(eigvects, 2), length(params.q), length(params.l))
    for i in 1:size(eigvects, 2)
        gipr[i, :, :] = compute_gipr(params, eigvects[:, i])
    end
    
    return compute_τ(params::MFAParameters, gipr::Array{F, 3})
end

function compute_ταf(params::MFAParameters, gipr::AbstractArray{F, 3}, μqlnμ::AbstractArray{F, 3}) where F
    gipr_mean = dropmean(gipr, dims=1)
    μqlnμ_mean = dropmean(μqlnμ, dims=1)
    τ = similar(gipr_mean)
    α = similar(gipr_mean)
    f = similar(gipr_mean)

    for i in 1:size(α, 2)
        ϵ = params.l[i]/params.ltc.N
        τ[:, i] = log.(gipr_mean[:, i]) ./ log(ϵ)
        α[:, i] = μqlnμ_mean[:, i] ./ (gipr_mean[:, i]*log(ϵ))
        f[:, i] = α[:, i].*params.q .- τ[:, i]
    end

    return τ, α, f
end

function compute_ταf(params::MFAParameters, eigvects::Array{F, 2}) where F
    gipr = Array{Float64}(undef, size(eigvects, 2), length(params.q), length(params.l))
    μqlnμ = similar(giprs)
    
    for i in 1:size(eigvects, 2)
        gipr[i, :, :], μqlnμ[i, :, :] = compute_gipr_2(params, @view eigvects[:, i])
    end
    
    return compute_ταf(params, gipr, μqlnμ)
end

end # module
