#Maximum enhancement

using JLD
using Plots
using Dates
using Polynomials
include("../library/abf2_pnscan.jl")
file = "/Users/pcs/data/ABF/rawdata/nu3-pn-break-t-sym/W01/pn.jld"

vars = load(file)

println(vars["info"])
W = collect(LinRange(vars["W_min"], vars["W_max"], vars["W_num"]))
ϕ = collect(LinRange(vars["phi_min"][1], vars["phi_max"][1], vars["phi_num"][1]))

PN_mean = vars["PN_mean"]
N = vars["N"]
PN_var = vars["PN_var"]
PN_hist = vars["PN_hist"]
plot(ϕ,PN_mean[:,1,15])
