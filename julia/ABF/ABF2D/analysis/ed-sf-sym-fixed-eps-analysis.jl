using DataFrames, CSV, MAT
using ROAG, Binning
using Statistics
using StatsBase
using Glob
using Plots
using LsqFit
using LaTeXStrings
##
len_W = 10
L = [40 50 60]
rdir = "/Users/pcs/codes/project/julia/ABF/ABF2D/V2test/"
# dir = ["L$(L[i])_Th1_W$(j)_E1.csv" for i in 1:length(L), j in 1:len_W]
dir = ["L$(L[i])_Th1_W1_E$(j).csv" for i in 1:length(L), j in 1:len_W]

ipr_mean = Array{Float64}(undef, len_W, length(L))
ipr_std = similar(ipr_mean)
ipr_ste = similar(ipr_mean)
for i in 1:length(L)
    for j in 1:len_W
        df = CSV.read(rdir*dir[i, j], DataFrame)
        ipr_mean[j, i] = mean(df.q2)
        ipr_std[j, i] = std(df.q2)
        ipr_ste[j, i] = ipr_std[j, i] / sqrt(length(df.E))
    end
end

plot(ipr_mean, yerror = ipr_ste)

τ = similar(ipr_mean)
τ_err = similar(ipr_mean)
for i in 1:length(L)
    τ[:, i] = log.(ipr_mean[:, i]) ./ log(0.05)
    τ_err[:, i] = abs.(ipr_ste[:, i] ./ ipr_mean[:,i] ./ log(0.05))
end

plot(τ, yerror = τ_err)
