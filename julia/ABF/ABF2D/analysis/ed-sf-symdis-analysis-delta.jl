using DataFrames, CSV, MAT
using ROAG, Binning
using Statistics
using StatsBase
using Glob
using Plots
using LsqFit
using LaTeXStrings
##
marker = (:circle, 4, 1., stroke(-0.5, 1., :black))
line = (:line, :solid, 1.5)
palette_roag = :Dark2_5
default(
    framestyle = :box,
    size = (600,400),
    # right_margin = [3mm 0mm],
    grid = false,
    minorticks = true,
    legend = (0.1, 0.75),
    fontfamily = "computer modern",
    tickfontsize = 13,
    guidefontsize = 13,
    legendfontsize = 13, palette = :default)

len_W = 1
len_E = 1
len_th = 10
L = [50 100 150]

rdir = ["/Users/pcs/data/ABF2D/symplectic/delta1_e0/L$l/1/" for l in L]
savedir = "/Users/pcs/data/ABF-sum/2d-sf-sym-pn-figs/"

dir = ["L$(L[i])_Th$(j)_W$(k)_E$(l).csv" for i in 1:length(L), j in 1:len_th, k in 1:len_W, l in 1:len_E]
ipr_mean = Array{Float64}(undef, length(L), len_th, len_W, len_E)
ipr_std = similar(ipr_mean)
ipr_ste = similar(ipr_mean)
E = similar(ipr_mean)
for i in 1:length(L), j in 1:len_th, k in 1:len_W, l in 1:len_E
    df = CSV.read(rdir[i]*dir[i, j, k, l], DataFrame)
    ipr_mean[i,j, k,l] = mean(df.q2)
    E[i, j, k, l] = mean(df.E)
    ipr_std[i, j, k, l] = std(df.q2)
    ipr_ste[i, j, k, l] = ipr_std[i, j, k, l] / sqrt(length(df.E))
end

plot(ipr_mean[:, :, 1, 1], yerror = ipr_ste[: ,:, 1, 1])

τ = similar(ipr_mean)
τ_err = similar(ipr_mean)
τ_std = similar(ipr_mean)

for i in 1:length(L), j in 1:len_th, k in 1:len_W, l in 1:len_E
    τ[i, j, k, l] = log.(ipr_mean[i, j, k, l]) ./ log(0.05)
    τ_err[i, j, k, l] = abs.(ipr_ste[i, j, k, l] ./ ipr_mean[i, j, k, l] ./ log(0.05))
    τ_std[i, j, k, l] = abs.(ipr_std[i, j, k, l] ./ ipr_mean[i, j, k, l] ./ log(0.05))
end

j = 2
p = plot(legend = :bottomright);
plot!(p, τ[:, :,  1, 1]', yerror = τ_err[:, :, 1, 1]')
xlabel!(L"E");
ylabel!(L"\tilde{\tau}")

# CSV.write(savedir*"clean_W0.01.csv", df)
# savefig(p, savedir*"clean_W0.01.pdf")
