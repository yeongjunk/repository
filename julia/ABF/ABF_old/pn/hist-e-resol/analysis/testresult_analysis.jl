#Maximum enhancement

using JLD
using Plots
using Dates

include("../library/abf2_pnscan.jl")
file = "/Users/pcs/Data/ABF/pn/rawdata/hist_e/025pi/pn.jld"

vars = load(file)

println(vars["info"])
W = collect(LinRange(vars["W_min"], vars["W_max"], vars["W_num"]))
θ = collect(LinRange(vars["theta_min"][1], vars["theta_max"][1], vars["theta_num"][1]))

PN_mean = vars["PN_mean"]
N = vars["N"]
PN_var = vars["PN_var"]
PN_hist = vars["PN_hist"]

j = 1
x = midpoints(PN_hist[j].edges[2])
y1 = PN_hist[j].weights[3,:]
y3 = PN_hist[j].weights[9,:]
y5 = PN_hist[j].weights[15,:]

E_bin_num = vars["E_bin_num"]
E_bins = range(-0.5, 0.5, length = E_bin_num)
E = (E_bins .+ 1) ./ W[j]

plot(x,(normalize(y1)),label="ΔE = $(round(E_bins[3],digits = 2))W")
plot!(x,(normalize(y3)), label="ΔE/W = $(round(E_bins[9],digits = 2))W")
plot!(x,normalize(y5), label="ΔE/W = $(round(E_bins[15],digits = 2))W")
vline!([4], color =:black, label = "CLS")
title!("Histogram")
xlabel!("PN")
ylabel!("Normalized Sample #")
plot(E_bins,PN_mean[1,1,:])
xlabel!("ΔE/W")
ylabel!("mean(PN)")
W[1]
