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
line = (:line, :solid, 2)
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

len_W = 20
L = [40 50 60 100]
rdir = "/Users/pcs/codes/project/julia/ABF/ABF2D/V2fine/"
savedir = "/Users/pcs/data/ABF-sum/2d-sf-sym-pn/"

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
τ_std = similar(ipr_mean)

for i in 1:length(L)
    τ[:, i] = log.(ipr_mean[:, i]) ./ log(0.05)
    τ_err[:, i] = abs.(ipr_ste[:, i] ./ ipr_mean[:,i] ./ log(0.05))
    τ_std[:, i] = abs.(ipr_std[:, i] ./ ipr_mean[:,i] ./ log(0.05))
end
lbl = [L"$L = %$(L[i])$" for i in 1:length(L)]
lbl = reshape(lbl, 1, length(L))

E= range(0.001, 0.17, length = 20)
p = plot(E, τ, yerror = τ_err, legend = :bottomleft, line = line, label = lbl)
xlabel!(L"E")
ylabel!(L"\tilde{\tau}")

p2 = plot(E, τ_std, legend = :bottomleft, line = line, label = lbl)


savefig(p, savedir*"tau_L40_100_th_0.25.pdf")
