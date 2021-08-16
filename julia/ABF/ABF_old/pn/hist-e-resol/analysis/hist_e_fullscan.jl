#Maximum enhancement

using JLD
using Plots
using Dates
using Polynomials
include("../library/abf2_pnscan.jl")
file = "/Users/pcs/Data/ABF/pn/rawdata/hist_e/fullscan/10/pn.jld"

vars = load(file)

println(vars["info"])
W = collect(LinRange(vars["W_min"], vars["W_max"], vars["W_num"]))
θ = collect(LinRange(vars["theta_min"][1], vars["theta_max"][1], vars["theta_num"][1]))

PN_mean = vars["PN_mean"]
N = vars["N"]
PN_var = vars["PN_var"]
PN_hist = vars["PN_hist"]

i = 5
j = 1
idx_E = [5 8 11 14 16]
x = midpoints(PN_hist[i,j].edges[2])
y = [PN_hist[i,j].weights[idx_E[k],:] for k in 1:length(idx_E)]

E_bin_num = vars["E_bin_num"]
E_bins = range(-0.5, 0.5, length = E_bin_num)
E = (E_bins .+ 1) ./ W[j]

p = plot()
for k = 1:length(idx_E)
    plot!(p, x,(normalize(y[k])),label="ΔE = $(round(E_bins[idx_E[k]],digits = 2))W")
end
display(p)
xlabel!("PN")
ylabel!("Probability Amplitude")

p_log = plot()
for k = 1:length(idx_E)
    plot!(p_log, x,log.(normalize(y[k])),label="ΔE = $(round(E_bins[idx_E[k]],digits = 2))W")
end
display(p_log)
xlabel!("PN")
ylabel!("ln(Probability Amplitude)")

idx_E_full = 8:24
y_full = ([log.(PN_hist[i,j].weights[idx_E_full[k],:]) for k in 1:length(idx_E_full)])
x[200]
x[320]

pfit = [Polynomials.fit(collect(x[200:340]), y_full[k][200:340], 1) for k in 1:length(idx_E_full)]
slope = Array{Float64}(undef, length(idx_E_full))
for k in 1:length(idx_E_full)
    slope[k] = pfit[k][1]
end

plot(E_bins[idx_E_full],slope)


p_loglog = plot()
1:length(idx_E)
    plot!(p_loglog, log.(x),log.(normalize(y[i])),label="ΔE = $(round(E_bins[idx_E[i]],digits = 2))W")
end
display(p_loglog)
xlims!(1, 3)


plot(x,log.(normalize(y1)),label="ΔE = $(round(E_bins[3],digits = 2))W")
plot!(x,log.(normalize(y3)), label="ΔE/W = $(round(E_bins[9],digits = 2))W")
plot!(x,log.(normalize(y5)), label="ΔE/W = $(round(E_bins[16],digits = 2))W")

plot(log.(x),log.(normalize(y1)),label="ΔE = $(round(E_bins[3],digits = 2))W")
plot!(log.(x),log.(normalize(y3)), label="ΔE/W = $(round(E_bins[9],digits = 2))W")
plot!(log.(x),log.(normalize(y5)), label="ΔE/W = $(round(E_bins[16],digits = 2))W")


# vline!([4], color =:black, label = "CLS")
title!("Histogram")
xlabel!("PN")
ylabel!("Normalized Sample #")
plot(E_bins,PN_mean[1,1,:])
xlabel!("ΔE/W")
ylabel!("mean(PN)")

θ

plotly()
