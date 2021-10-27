using DataFrames, CSV, MAT
using ROAG, Binning
using Statistics
using StatsBase
using Glob
using Plots
using LsqFit
using LaTeXStrings


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

##
L = [100 200 300 400 500]
rdir = "/Users/pcs/data/ABF2D/pd-fixed-eps/"
dir = [rdir*"L$(L[i])/" for i in 1:length(L)]
#
# k = 6
# for i in 2:7
#     init = (i-1)*3
#     for j in 1:3
#         mv(dir[k]*"$(i)/L$(L[k])_Th1_W1_E$(j).csv", dir[k]*"$(i)/L$(L[k])_Th1_W1_E$(init + j).csv")
#     end
# end

E_c = range(0., 0.3, length = 21)
θ = 0.25

savedir = "/Users/pcs/data/ABF-sum/2d-sf-pd-pn/"

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
p = plot()
for i in 1:length(L)
    p = plot!(E_c, τ[:, i], yerror = τ_err, label = L"L = %$(L[i])", legend = :bottomright, palette = :tab10, line = line)
end
xlabel!(L"E")
ylabel!(L"\tilde{\tau}")
savefig(p, savedir*"tau_th025.pdf")
xlims!(-0.01, 0.05)
ylims!(1.89, 1.92)
