using LinearAlgebra, SparseArrays
using Random
using DataFrames, CSV
using ArgParse, JSON

include("./chain-ed-arpack.jl") # read parameters from configuration file

function readconfig(d::Dict)
    R = d["R"]
    J = d["J"]
    N = d["N"]
    nev = d["nev"]
    ga = d["ga"]
    W = d["W"]
    λ2 = collect(d["lam_min"]:d["lam_step"]:d["lam_max"]).^2
    return ScanParameters(R = R, J = J, N = N, nev = nev, ga = ga, W = W, λ2 = λ2)
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
    fn = "bdgchain-W$(p.W)-ga$(p.ga)-nev$(p.nev)-R$(p.R).csv"
    @time λ, pn = scan_bdg(p)

    CSV.write(fn, DataFrame(λ = λ, pn = pn))
end

main(ARGS)
