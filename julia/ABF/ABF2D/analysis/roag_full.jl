using DataFrames, CSV
modulepath =  "/Users/pcs/codes/chain/ABF/module"
push!(LOAD_PATH, modulepath)
using ROAG, Binning
using Statistics
using StatsBase
using MAT
using Plots, LaTeXStrings
using Plots.PlotMeasures
## Function definitions
function load_file(dir, i, j, L, R, θ)
    subdir = "L$(L[i])/"
    fn = "L$(L[i])_Th$(j)_R$(R[i]).csv"
    return Tables.columntable(CSV.read(dir*subdir*fn, DataFrame))
end

function roag_mean_binned(E::AbstractArray{T}, E_edges::AbstractArray{T}) where T <: Number
    r_mean = Vector{Float64}(undef, length(E_edges)-1)
    lbl = binning_unsorted(E, E_edges) # binning
    for i in 1:(length(E_edges) - 1)
        r = E[findall(a -> a == i, lbl)]
        if length(r) < 3
            r_mean[i] = NaN
        else
            roag!(r)
            r_mean[i] = mean(r)
        end
    end
    return r_mean::Vector{Float64}
end

function roag_scan(;dir, L, R, θ, E_edges)
    r_m_m = Array{Float64}(undef, length(E_edges)-1 ,length(θ), length(L))
    r_m_std = similar(r_m_m)
    for i in 1:length(L)
        println("i = ", i)
        for j in 1:length(θ)
            df = load_file(dir, i, j, L, R, θ)
            r_m = Array{Float64}(undef, length(E_edges)-1, R[i])
            for r in 1:R[i]
                idx = searchsorted(df.r, r)
                E = df.E[idx]
                r_m[:, r] .= roag_mean_binned(E, E_edges)
            end
            for k in 1:length(E_edges) - 1
                r_m_m[k, j, i] = mean(filter(!isnan, r_m[k, :]))
                r_m_std[k, j, i] = std(filter(!isnan, r_m[k, :]))
            end
        end
    end
    return r_m_m, r_m_std
end

## Parameters
L = [20 40 60 80 100]
R = [10000 2500 1000 1000 500]
θ = range(0.0001, 0.25, length = 26)
E_edges = collect(0.65:0.1:1.35)
pushfirst!(E_edges,0.60)
push!(E_edges,  1.4)

##
dir = "/Users/pcs/data/ABF2D/onsite/"
savedir = "/Users/pcs/data/ABF2D/onsite/analyzed-data/"
##
r_m_m, r_m_std = roag_scan(dir = dir, L = L, R = R , θ = θ, E_edges = E_edges)

matwrite(savedir*"roag.mat",
    Dict("r_m_m" => r_m_m,
    "r_m_std" => r_m_std)) # note that this does NOT introduce a variable ``varname`` into scope
