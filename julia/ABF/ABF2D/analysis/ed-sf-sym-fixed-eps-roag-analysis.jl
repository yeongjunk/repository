using DataFrames, CSV, MAT
using ROAG, Binning
using Statistics
using StatsBase
using Glob
using Plots
using LsqFit
using LaTeXStrings
using ROAG
##
# marker = (:circle, 4, 1., stroke(0, 1., :black))
# line = (:line, :solid, 1.5)
# palette_roag = :Dark2_5
# default(
#     framestyle = :box,
#     size = (600, 400),
#     # right_margin = [3mm 0mm],
#     grid = false,
#     minorticks = true,
#     # legend = (0.1, 0.75),
#     fontfamily = "computer modern",
#     tickfontsize = 13,
#     guidefontsize = 13,
#     legendfontsize = 13, palette = :default)

len_W = 20
L = [50 100 200]
rdir = ["/Users/pcs/data/ABF2D/symplectic/L$(l)/" for l in L]
savedir = "/Users/pcs/data/ABF-sum/2d-sf-sym-pn/"

# dir = ["L$(L[i])_Th1_W$(j)_E1.csv" for i in 1:length(L), j in 1:len_W]
dir = ["L$(L[i])_Th4_W1_E$(j).csv" for i in 1:length(L), j in 1:len_W]

roag_mean = Array{Float64}(undef, len_W, length(L))
roag_std = similar(roag_mean)
roag_ste = similar(roag_mean)
E_band = similar(roag_mean)
for i in 1:length(L)
    for j in 1:len_W
        df = CSV.read(rdir[i]*dir[i, j], DataFrame)
        roag_means = Float64[]
        E_band[j, i] = df.E[1]
        for R in 1:maximum(df.r)
            df_r = df[df.r .== R, :]
            E = copy(df_r.E)
            sort!(E)
            E = E[1:2:end]
            roag!(E)
            r_mean = mean(E)
            push!(roag_means, mean(E))
            if isnan(r_mean)
                error("There is NaN i = $(i), j = $(j)")
            end
        end
        roag_mean[j, i] = mean(roag_means)
        roag_std[j, i] = std(roag_means)
        roag_ste[j, i] = roag_std[j, i] / sqrt(length(roag_means))
    end
end
plotly()

p = plot()
for i in 1:length(L)
    plot!(p, E_band[:, i], roag_mean[:, i], yerror = roag_ste)
end
xlabel!("E")
ylabel!("r")

savefig(p, savedir*"roag.pdf")
