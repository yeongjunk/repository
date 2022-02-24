using CSV
using DataFrames
using LinearAlgebra
using MAT
using JSON
using Plots
using Glob
using Statistics
using MAT

global opendir = "/Users/pcs/data/ABF/rawdata/phase/1d-sf-pd-pn/"
global root_dir = @__DIR__
include(root_dir*"/params.jl")
global L = [300, 1000, 1001, 5000]
global dir = [opendir for i in 1:length(L)]
config_lists = opendir*"config_L5000"
config  = JSON.parsefile(config_lists)

global th = range(config["th"][1], config["th"][2], length = config["th"][3])
global W = range(config["W"][1], config["W"][2], length = config["W"][3])
global E = range(config["E"][1], config["E"][2], length = config["E"][3])
global V = config["V"]

function opener(i, j, k, l)
    #----------- Read configuration files ---------#
    #---------- open data ---------#
    fn = "L$(L[i])_Th$(j)_W$(lpad(k, 2, "0"))_E$(lpad(l, 2, "0")).csv"
    df = CSV.read(dir[i]*fn, DataFrame)
    return df, L[i], th[j], W[k], E[l]
end

function processor_ipr(b; q = 2)
    ipr_mean = Array{Float64}(undef, length(L), length(th), length(W), length(E))
    ipr_std = similar(ipr_mean)
    ipr_ste = similar(ipr_mean)
    e = similar(ipr_mean)
    for i in 1:length(L), j in 1:length(th), k in 1:length(W), l in 1:length(E)
        df, _, _, _  = opener(i, j, k, l, printinfo = false)
        ipr_mean[i,j, k,l] = mean(df[:,"l$(b[i])_q$(q)"])
        e[i, j, k, l] = mean(df.E)
        ipr_std[i, j, k, l] = std(df[:, "l$(b[i])_q$(q)"])
        ipr_ste[i, j, k, l] = ipr_std[i, j, k, l] / sqrt(length(df.E))
    end
    return e, ipr_mean, ipr_std, ipr_ste
end

function processor_tau(b; q = 2)
    e, ipr_mean, ipr_std, ipr_ste = processor_ipr(b, q=q)

    τ = similar(ipr_mean)
    τ_err = similar(ipr_mean)
    τ_std = similar(ipr_mean)

    for i in 1:length(L), j in 1:length(th), k in 1:length(W), l in 1:length(E)
        τ[i, j, k, l] = log(ipr_mean[i, j, k, l]) / log(b[i]/L[i])
        τ_err[i, j, k, l] = abs(ipr_ste[i, j, k, l] / ipr_mean[i, j, k, l] / log(b[i]/L[i]))
        τ_std[i, j, k, l] = abs(ipr_std[i, j, k, l] / ipr_mean[i, j, k, l] / log(b[i]/L[i]))
    end
    return e, τ, τ_err, τ_std
end

l = 1
q = 2
ipr = Array{Float64}(undef, length(L), length(th), length(E))
ipr_std = similar(ipr)
ipr_ste = similar(ipr)
ipr_E = similar(ipr)
for i in 1:length(L), j in 1:length(th), k in 1:length(E)
    df, _, _, _ = opener(i, j, 1, k)
    ipr_E[i, j, k] = mean(df.E)
    ipr[i, j, k] = mean(df[:, "l$(l)_q$(q)"])
    ipr_std[i, j, k] = std(df[:, "l$(l)_q$(q)"])
    ipr_ste[i, j, k] = ipr_std[i, j, k]./length(df[:, "l$(l)_q$(q)"])
end
plot(ipr_E[:, 2, :]', (ipr[:, 2, :]').^-1)

L = L
th = vec(collect(th))
W = vec(collect(W))
V = V
description = "IPR of Scale free model of d = 1 nu = 2 ABF. First axis: L, Second axis: thta, Thrid axis: E"

data = Dict("L" => L,
    "th" => th, "V" => V, "description" => description,
    "E" => ipr_E, "ipr" => ipr, "ipr_std" => ipr_std,
    "ipr_ste" => ipr_ste)

matwrite(root_dir*"/1d-sf-phase-ipr.mat", data, compress=true)
