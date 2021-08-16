#Maximum enhancement

using JLD
using Plots
using Dates
using StatsBase

include("../library/abf2_pnscan.jl")
dir = "/Users/pcs/data/ABF/rawdata/hop_dis_weak"

vars = load(dir*"/pn.jld")

println(vars["info"])

p = readconfig(vars)
V, θ, _ = expand_params(p)

PN_mean = vars["PN_mean"]
N = vars["num_samp"]
PN_var = vars["PN_var"]
PN_hist = vars["PN_hist"]


#PN_mean
p1 = plot()
for i in 1:length(V)
    plot!(θ,PN_mean[:,i,10], label = "delta = $(round(10 .^V[i]/p.W, digits = 4))", legend = :topleft)
end
display(p1)
xlabel!("θ/π")
ylabel!("<PN>")
savefig(p1,"/Users/pcs/data/ABF/analysis/nu2-pn-phasedis/mean.png")
p2 = plot()
for i in 1:length(V)
    yy= normalize(PN_hist[25,i].weights[10,:])
    xx = midpoints(PN_hist[25,i].edges[2])
    plot!(xx,log10.(yy), label = "δ = $(round(10 .^V[i]/p.W, digits = 4))")
end
xlabel!("PN")
ylabel!("Probability")
savefig(p2,"/Users/pcs/data/ABF/analysis/nu2-pn-phasedis/histlog.png")

display(p2)
