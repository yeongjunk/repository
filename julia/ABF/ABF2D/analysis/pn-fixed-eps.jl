using DataFrames, CSV, MAT
using ROAG, Binning
using Statistics
using StatsBase
using Glob
using Plots
using LsqFit
using LaTeXStrings



##
L = [100 150 200 250 300 400 500]
rdir = "/Users/pcs/data/ABF2D/pd-fixed-eps/W01/"
dir = [rdir*"L$(L[i])/" for i in 1:length(L)]

# k = 7
# for i in 2:7
#     init = (i-1)*3
#     for j in 1:3
#         mv(dir[k]*"$(i)/L$(L[k])_Th1_W1_E$(j).csv", dir[k]*"$(i)/L$(L[k])_Th1_W1_E$(init + j).csv")
#     end
# end

E_c = range(0., 0.3, length = 21)
θ = 0.25

savedir = "/Users/pcs/data/ABF2D/analyzed-data/"

ipr_mean = Array{Float64}(undef, length(E_c), length(L))
ipr_std = similar(ipr_mean)
ipr_ste = similar(ipr_mean)
for i in 1:length(L)
    for j in 1:length(E_c)
        df = CSV.read(dir[i]*"L$(L[i])_Th1_W1_E$(j).csv", DataFrame)
        ipr_mean[j, i] = mean(df.q2)
        ipr_std[j, i] = std(df.q2)
        ipr_ste[j, i] = ipr_std[j, i] / sqrt(length(df.E))
    end
end

τ = similar(ipr_mean)
τ_err = similar(ipr_mean)
for i in 1:length(L)
    τ[:, i] = log.(ipr_mean[:, i]) ./ log(0.1)
    τ_err[:, i] = abs.(ipr_ste[:, i] ./ ipr_mean[:,i] ./ log(0.1))
end

plot(E_c, τ, yerror = τ_err, label = "L = ".*string.(L), legend = :bottomright, palette = :tab10)

xlims!(-0.01, 0.05)
ylims!(1.89, 1.92)
