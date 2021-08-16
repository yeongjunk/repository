using DataFrames, CSV
using Statistics, StatsBase
using Plots, Interpolations, LaTeXStrings
using ROAG, Binning # Custom modules


##
L = [100 1000]
R = [10000 1000]
θ = range(0.0001, 0.25, length = 26)
E_edges = collect(0.65:0.1:1.35)
pushfirst!(E_edges,0.60)
push!(E_edges,  1.4)

##
dir = "/Users/pcs/data/ABF1D/onsite/"
savedir = "/Users/pcs/data/ABF1D/onsite/analyzed-data/"
#
##
function load_file(dir, i, j, L, R, θ)
    subdir = "L$(L[i])/"
    fn = "L$(L[i])_Th$(j)_R$(R[i]).csv"
    return Tables.columntable(CSV.read(dir*subdir*fn, DataFrame))
end

function pn_mean_binned(E::AbstractArray{T}, pn::AbstractArray{T}, E_edges::AbstractArray{T}) where T <: Number
    pn_mean = Vector{Float64}(undef, length(E_edges)-1)
    lbl = binning_sorted(E, E_edges) # binning
    for i in 1:length(E_edges) - 1
        @views pn_mean[i] = mean(pn[findall(a -> a == i, lbl)])
    end
    return pn_mean::Vector{Float64}
end

function pn_scan(;dir, L, R, θ, E_edges)
    pn_m_m = Array{Float64}(undef, length(E_edges)-1 ,length(θ), length(L))
    pn_m_std = similar(pn_m_m)
    for i in 1:length(L)
        println("i = ", i)
        for j in 1:length(θ)
            df = load_file(dir, i, j, L, R, θ)
            pn_m = Array{Float64}(undef, length(E_edges)-1, R[i])
            for r in 1:R[i]
                idx = searchsorted(df.r, r)
                @views E = df.E[idx]
                @views pn = (df[3][idx]).^-1
                pn_m[:, r] .= pn_mean_binned(E, pn, E_edges)
            end
            for k in 1:length(E_edges) - 1
                pn_m_m[k, j, i] = mean(filter(!isnan, pn_m[k, :]))
                pn_m_std[k, j, i] = std(filter(!isnan, pn_m[k, :]))
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
    Y = log.(pn_m_m[i,j,:])
    v = X\Y
    pn_exp_fit[i,j] = v[2]
end

pn_exp_diff = Array{Float64}(undef, length(E_edges)-1, length(θ), length(L)-1)
for i = 1:length(E_edges)-1, j = 1:length(θ)
    dx = diff(log.(vec(L)))
    dy = diff(log.(pn_m_m[i,j,:]))
    pn_exp_diff[i,j,:] = dy./dx
end

matwrite(savedir*"pn.mat",
    Dict("pn_m_m" => pn_m_m,
    "pn_m_std" => pn_m_std,
    "pn_exp_fit" => pn_exp_fit,
    "pn_exp_diff" => pn_exp_diff))

    itp = interpolate(pn_m_m[:,:,end], BSpline(Linear()))
    r_new = itp(range(1, length(E_edges)-1, length = 100), range(1,26, length = 100))
    θ_itp = interpolate(θ, BSpline(Linear()))
    θ_new  = θ_itp(range(1,length(θ), length = 100))
    E_itp = interpolate(midpoints(E_edges), BSpline(Linear()))
    E_new = E_itp(range(1,length(E_edges)-1, length = 100))


marker = (:circle, 4, 1., stroke(-0.1, 1., :black))
line = (:line, :solid, 2)
palette_roag = :default
default(
    framestyle = :box,
    size = (600,400),
    # right_margin = [3mm 0mm],
    grid = false,
    minorticks = true,
    legend = (0.1, 0.75),
    fontfamily = "computer modern",
    tickfontsize = 14,
    guidefontsize = 15,
    legendfontsize = 10)

p = plot(
    xlabel = L"$\Theta / \pi$",
    ylabel = L"$PN$")
for i in 1:2
    plot!(θ, pn_m_m[5, :, i], line = line, marker = marker, label = L"$L = %$(L[i])$")
end
annotate!(0.02, 90, L"E = 1")
p
