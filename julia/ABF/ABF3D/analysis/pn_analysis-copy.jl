using DataFrames, CSV, MAT
using ROAG, Binning
using Statistics
using StatsBase


##
L = [10 15 20 25 30]
R = [1200 400 100 65 50]
dir = ["/Users/pcs/data/ABF-sum/3d-sf-on-pn-roag/full-ipr/L$(L[i])/" for i in 1:length(L)]
θ = range(0.0001, 0.25, length = 26)
savedir = "/Users/pcs/data/ABF3D/analyzed-data/"
E_edges = collect(range(0.6, 1.4, length = 40))
# E_edges = collect(0.65:0.1:1.35)
# pushfirst!(E_edges,0.625)
# pushfirst!(E_edges, 0.6)
# push!(E_edges,  1.375)
# push!(E_edges,  1.4)
# push!(E_edges, 1.30)
# push!(E_edges, 0.70)
# sort!(E_edges)
##
function load_file(dir, i, j, L, R, θ)
    fn = "L$(L[i])_Th$(j)_R$(R[i]).csv"
    return CSV.read(dir*fn, DataFrame)
end

function pn_mean_binned(E::AbstractArray{T}, pn::AbstractArray{T}, E_edges::AbstractArray{T}) where T <: Number
    pn_mean = Vector{Float64}(undef, length(E_edges)-1)
    lbl = binning_sorted(E, E_edges) # binning
    for i in 1:length(E_edges) - 1
        @views pn_mean[i] = mean(log.(pn[findall(a -> a == i, lbl)]))
    end
    return pn_mean::Vector{Float64}
end

function pn_scan(;dir, L, R, θ, E_edges)
    pn_m_m = Array{Float64}(undef, length(E_edges)-1 ,length(θ), length(L))
    pn_m_std = similar(pn_m_m)
    for i in 1:length(L)
        println("i = ", i)
        for j in 1:length(θ)
            df = load_file(dir[i], i, j, L, R, θ)
            pn_m = Array{Float64}(undef, length(E_edges)-1, R[i])

            for r in 1:R[i]
                idx = searchsorted(df.r, r)
                E = df.E[idx]
                pn = 1. ./ (df[idx, Symbol("2.0")])
                pn_m[:, r] .= pn_mean_binned(E, pn, E_edges)
            end
            for k in 1:length(E_edges) - 1
                pn_m_m[k, j, i] = mean(pn_m[k, :])
                pn_m_std[k, j, i] = std(pn_m[k, :])
            end
        end
    end
    return pn_m_m, pn_m_std
end
# 0.6 < E < 1.4 (slightly smaller than (1-W/2) and (1 + W/2). The midband is E = 1.
## ROAG for the midband: scaling of ROAG.
pn_m_m, pn_m_std = pn_scan(dir = dir, L = L, R = R , θ = θ, E_edges = E_edges)

pn_exp_fit = Array{Float64}(undef, length(E_edges)-1, length(θ))
for i = 1:length(E_edges)-1, j = 1:length(θ)
    X = [ones(length(L)) log.(vec(L))]
    Y = pn_m_m[i,j,:]
    v = X\Y
    pn_exp_fit[i,j] = v[2]
end

pn_exp_diff = Array{Float64}(undef, length(E_edges)-1, length(θ), length(L)-1)
for i = 1:length(E_edges)-1, j = 1:length(θ)
    dx = diff(log.(vec(L)))
    dy = diff(pn_m_m[i,j,:])
    pn_exp_diff[i,j,:] = dy./dx
end

matwrite(savedir*"pn_log.mat",
    Dict("pn_m_m" => pn_m_m,
    "pn_m_std" => pn_m_std,
    "pn_exp_fit" => pn_exp_fit,
    "pn_exp_diff" => pn_exp_diff))
