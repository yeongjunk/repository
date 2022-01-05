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
    tickfontsize = 15,
    guidefontsize = 20,
    legendfontsize = 15, palette = :default)

len_W = 10
L = [50 100 200 300]
# rdir = ["/Users/pcs/data/ABF-sum/2d-sf-sym-pn-2/L$(l)/" for l in L]
# savedir = "/Users/pcs/data/ABF-sum/2d-sf-sym-pn-figs/"
rdir = ["/Users/pcs/codes/project/julia/ABF/ABF2D/" for l in L]
# rdir = ["/Users/pcs/data/ABF-sum/2d-fe-clean-sym/" for l in L]
savedir = "/Users/pcs/data/ABF-sum/2d-sf-sym-pn-figs/"

dir = ["L$(L[i])_Th4_W$(j)_E1.csv" for i in 1:length(L), j in 1:len_W]
# dir = ["L$(L[i])_Th4_W1_E$(j).csv" for i in 1:length(L), j in 1:len_W]

ipr_mean = Array{Float64}(undef, len_W, length(L))
ipr_std = similar(ipr_mean)
ipr_ste = similar(ipr_mean)
E = similar(ipr_mean)
for i in 1:length(L)
    for j in 1:len_W
        df = CSV.read(rdir[i]*dir[i, j], DataFrame)
        ipr_mean[j, i] = mean(df.q2)
        E[j, i] = mean(df.E)
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

p = plot(legend = :bottomleft);
for i in 1:length(L)
    plot!(p, E[:, i], τ[:, i], yerror = τ_err, line = line, label = L"$L = %$(L[i])$")
end
xlabel!(L"E");
ylabel!(L"\tilde{\tau}")

df = DataFrame(E = E[:, 1], tau = τ[:, 1])
CSV.write(savedir*"clean_W0.01.csv", df)
savefig(p, savedir*"clean_W0.01.pdf")
