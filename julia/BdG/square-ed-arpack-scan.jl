using LinearAlgebra, SparseArrays
using Random
using DataFrames, CSV
using ArgParse, JSON

include("./square-ed-arpack.jl") # read parameters from configuration file

function readconfig(d::Dict)
    R = d["R"]
    J = d["J"]
    N = d["N"]
    nev = d["nev"]
    ga = d["ga"]
    W = d["W"]
    λ2 = collect(d["lam_min"]:d["lam_step"]:d["lam_max"]).^2
    return [ScanParameters(R = R, J = J, N = N[i], nev = nev, ga = ga, W = W, λ2 = λ2) for i in 1:length(N)]
end

function main(ARGS)
    opts = ArgParseSettings(description="Scan and compute pn for all parameters of nu=2 ABF")
    @add_arg_table! opts begin
    "c"
        help = "configuration"
        arg_type = AbstractString
    end
    # Parse the arguments
    popts   = parse_args(opts)
    config  = JSON.parsefile(popts["c"])
    p = readconfig(config)
    for i in 1:length(p)
        fn = "bdgsquare-W$(p[i].W)-ga$(p[i].ga)-nev$(p[i].nev)-R$(p[i].R)-L$(p[i].N).csv"
        @time λ, pn = scan_bdg(p[i])
        CSV.write(fn, DataFrame(λ = λ, pn = pn))
    end
end

main(ARGS)
