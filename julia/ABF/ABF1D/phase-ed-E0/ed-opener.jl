using CSV
using DataFrames
using LinearAlgebra
using MAT
using JSON
using Plots
using Glob
using Statistics
using Binning

global opendir = "/Users/pcs/data/ABF-sum/raw-data/phase/1d-sf-pd-pn/"
global root_dir = @__DIR__
include(root_dir*"/lanczos-params.jl")
global L = vec([300 1000 5000])
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
ipr = Array{Float64}(undef, 3, 31)
E = similar(ipr)
for i in 1:3, j in 3:31
    df, _, _, _ = opener(i, 2, 1, j-1)
    E[i, j] = mean(df.E)
    ipr[i,j] = mean(df[:, "l$(l)_q$(q)"])
end

for i in 1:3, j in 1
    df, _, _, _ = opener(i, 2, 1, j)
    E_max = maximum(df.E)
    E_edges = vec([-E_max/100, E_max/100, E_max])
    labels = binning_unsorted(df.E, E_edges)
    idx_1 = findall(x -> x == 1, labels)
    idx_2 = findall(x -> x == 2, labels)
    E[i, 1] = mean(df.E[idx_1])
    E[i, 2] = mean(df.E[idx_2])
    ipr[i,1] = mean(df[idx_1, "l$(l)_q$(q)"])
    ipr[i,2] = mean(df[idx_2, "l$(l)_q$(q)"])
end

scatter(E, (ipr.^-1), legend = false)



df, _, _, _ = opener(1, 2, 1, 1)
E_max = maximum(df.E)
E_edges = vec([-E_max/1000, E_max/1000, E_max/10, E_max])
labels = binning_unsorted(df.E, E_edges)
idx_1 = findall(x -> x == 1, labels)
idx_2 = findall(x -> x == 2, labels)
idx_3 = findall(x -> x == 3, labels)
mean(df[idx_1,"l1_q2"])^-1
mean(df[idx_2,"l1_q2"])^-1
mean(df[idx_3,"l1_q2"])^-1


scatter(df.E[1:1000], df[1:1000, "l1_q2"].^ -1)
