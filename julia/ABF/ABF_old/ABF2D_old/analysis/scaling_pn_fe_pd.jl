using DataFrames
using Plots
using LinearAlgebra
using StatsBase
using CSV
const modulepath =  "/Users/pcs/codes/chain/ABF2D/module"
push!(LOAD_PATH, modulepath)
using ABF2D
using Lattice
using PN
using Statistics


raw_dir = "/Users/pcs/data/ABF2D/pd-fe2/"
sizes = ["L20" "L40" "L60" "L80" "L100"]
sizedir = raw_dir.*sizes
sizes_int = collect(20:20:100)
function readdir_csv(dir)
    files = readdir(dir, join = true)
    filter!(s->occursin(r".csv", s), files)
end

function cut_mean(E, PN, r)
    idx = findall(x -> abs(x) < r*E[end], E)
    return mean(PN[idx])
end

## Read all filenames
fns = Dict()
for i in 1:length(sizes)
    fns[sizes[i]] = readdir_csv(sizedir[i])
end

## Process the data
r = 0.05
PN_mean = Dict()
PN_mean_mean = Float64[]

for i in 1:length(sizes)
    PN_mean[sizes_int[i]] = Float64[]
    for j in 1:length(fns[sizes[i]])
        df = CSV.read(fns[sizes[i]][j], DataFrame)
        df.E .-= 1.
        push!(PN_mean[sizes_int[i]], cut_mean(df.E, df.PN, r))
    end
    push!(PN_mean_mean, mean(PN_mean[sizes_int[i]]))
end
display(PN_mean_mean)

## Linear Fitting
X = hcat(ones(length(sizes)), log10.(sizes_int))
Y = log10.(PN_mean_mean)

v = X\Y

## Plotting
scatter(log10.(sizes_int), log10.(PN_mean_mean))
plot!(log10.(sizes_int), v[1] .+ v[2]*log10.(sizes_int))

println(v[2])
