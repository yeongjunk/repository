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

len_W = 10
L = [50 100 150]
# rdir = ["/Users/pcs/data/ABF-sum/2d-sf-sym-pn-2/L$(l)/" for l in L]
# savedir = "/Users/pcs/data/ABF-sum/2d-sf-sym-pn-figs/"
rdir = ["/Users/pcs/data/ABF2D/symplectic/delta1/L$l/1/" for l in L]
# rdir = ["/Users/pcs/data/ABF-sum/2d-fe-clean-sym/" for l in L]
savedir = "/Users/pcs/data/ABF-sum/2d-sf-sym-pn-figs/"

#dir = ["L$(L[i])_Th4_W$(j)_E1.csv" for i in 1:length(L), j in 1:len_W]
θ = 1:4
dir = ["L$(L[i])_Th$(k)_W1_E$(j).csv" for i in 1:length(L),k in 1:length(θ), j in 1:len_W]
ipr_mean = Array{Float64}(undef, len_W, length(L), length(θ))
ipr_std = similar(ipr_mean)
ipr_ste = similar(ipr_mean)
E = similar(ipr_mean)
for i in 1:length(L)
    for k in 1:length(θ)
        for j in 1:len_W
            df = CSV.read(rdir[i]*dir[i, k, j], DataFrame)
            ipr_mean[j, i, k] = mean(df.q2)
            E[j, i, k] = mean(df.E)
            ipr_std[j, i, k] = std(df.q2)
            ipr_ste[j, i, k] = ipr_std[j, i, k] / sqrt(length(df.E))
        end
    end
end

plot(ipr_mean[:,:, 2], yerror = ipr_ste[:, :, 2])

τ = similar(ipr_mean)
τ_err = similar(ipr_mean)
τ_std = similar(ipr_mean)

for i in 1:length(L)
    for j in 1:length(θ)
    τ[:, i, j] = log.(ipr_mean[:, i, j]) ./ log(0.05)
    τ_err[:, i, j] = abs.(ipr_ste[:, i, j] ./ ipr_mean[:,i, j] ./ log(0.05))
    τ_std[:, i, j] = abs.(ipr_std[:, i, j] ./ ipr_mean[:,i, j] ./ log(0.05))
    end
end

j = 2
p = plot(legend = :bottomleft);
for i in 1:length(L)
    plot!(p, E[:, i, j], τ[:, i, j], yerror = τ_err[:, :, j], line = line, label = L"$L = %$(L[i])$")
end
xlabel!(L"E");
ylabel!(L"\tilde{\tau}")

df = DataFrame(E = E[:, 1], tau = τ[:, 1])
# CSV.write(savedir*"clean_W0.01.csv", df)
# savefig(p, savedir*"clean_W0.01.pdf")
