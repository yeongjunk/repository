#Maximum enhancement

using JLD
using Plots
using Dates
include("../library/abf2_pnscan.jl")
dir = "/Users/pcs/data/ABF/rawdata/nu2-pn-hist-e-resolved/fullscan/"

#Empty array
N = [];PN_var = [];PN_mean = [];PN_hist = []; θ = [];


for i in 1:10
    file = dir*"$(i)/pn.jld"
    vars = JLD.load(file)
    tPN_mean = vars["PN_mean"]
    tN = vars["N"]
    tPN_var = vars["PN_var"]
    tPN_hist = vars["PN_hist"]
    tθ = collect(LinRange(vars["theta_min"][1], vars["theta_max"][1], vars["theta_num"][1]))
    push!(N, tN)
    push!(PN_mean, tPN_mean)
    push!(PN_var, tPN_var)
    push!(PN_hist, tPN_hist)
    push!(θ, tθ)
end

N = reduce(vcat,N)
PN_mean = reduce(vcat,PN_mean)
PN_var = reduce(vcat,PN_var)
PN_hist = reduce(vcat,PN_hist)
θ = reduce(vcat,θ)


PN_cls = 4*(1 .+ cospi.(2θ).^2).^(-2)
PN_cls = hcat(θ,PN_cls)
PN_weak_pyform = Array{Float64}(undef,50, 2)

PN_weak_pyform[:,1] = θ
PN_weak_pyform[:,2] = PN_mean[:,1,16]

PN_weak_var_pyform = Array{Float64}(undef,50, 2)

PN_weak_var_pyform[:,1] = θ
PN_weak_var_pyform[:,2] = PN_var[:,1,16]

PN_hist_weak = PN_hist[:,1]

binlen_PN = length(PN_hist_weak[50].weights[16,:])
PN_weak_hist_pyform = Array{Float64}(undef,binlen_PN, 2)

PN_weak_hist_pyform[:,1] = collect(midpoints(PN_hist_weak[50].edges[2]))
PN_weak_hist_pyform[:,2] = normalize(PN_hist_weak[50].weights[16,:])


matdata = Dict("PN_mean" => PN_weak_pyform,
"PN_var" => PN_weak_var_pyform,
"PN_hist" => PN_weak_hist_pyform,
"PN_cls" => PN_cls
)

fn_out = dir*"/pn_weak_py.mat"
using MAT
matwrite(fn_out, matdata, compress=true)
println("saved")
