#Maximum enhancement

using JLD
using Plots
using Dates
using StatsBase

include("../library/abf2_pnscan.jl")

topdir = "/Users/pcs/data/ABF/rawdata/nu2-pn/phase-dis-eff-absbw/"
subdir = string.(collect(1:14))

#Empty data lists
θ = []
δ = similar(θ)
PN_mean = similar(θ)
PN_var = similar(θ)
PN_hist = similar(θ)
N = similar(θ)
E_bw = similar(θ)


δ_min = 0.
δ_max = 0.
δ_num = 0


for i in 1:length(subdir)
    dir = topdir*subdir[i]
    vars = load(dir*"/pn.jld")
    p = readconfig(vars)

    if i == 1
        global δ_min += p.δ_min
    end
    if i == length(subdir)
        global δ_max += p.δ_max
    end
    global δ_num += p.δ_num

    tθ, tδ = expand_params(p)
    tPN_mean = vars["PN_mean"]
    tN = vars["num_samp"]
    tPN_var = vars["PN_var"]
    tPN_hist = vars["PN_hist"]
    tE_bw = vars["E_bw"]

    push!(δ, tδ)
    push!(PN_mean,tPN_mean)
    push!(PN_var,tPN_var)
    push!(PN_hist,tPN_hist)
    push!(N,tN)
    push!(E_bw, tE_bw)
end

function deltacat(A,B)
    cat(A,B, dims = 2) #concatnate along second dim
end

N = reduce(deltacat,N)
PN_mean = reduce(deltacat,PN_mean)
PN_var = reduce(deltacat,PN_var)
PN_hist = reduce(deltacat,PN_hist)
E_bw = reduce(deltacat, E_bw)
δ = reduce(vcat,δ)


dir = topdir*subdir[1]
vars = load(dir*"/pn.jld")
vars["delta_min"] = δ_min
vars["delta_max"] = δ_max
vars["delta_num"] = δ_num

vars["PN_mean"] = PN_mean
vars["num_samp"] = N
vars["PN_var"] = PN_var
vars["PN_hist"] = PN_hist
vars["E_bw"] = E_bw

fn_out = "pn_merged.jld"
save(topdir*fn_out, vars)
print("Merged")
