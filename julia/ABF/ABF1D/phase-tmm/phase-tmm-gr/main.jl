using Distributed
using JSON, ArgParse
using CSV
using DataFrames
using SharedArrays
using DelimitedFiles
@everywhere include("./xi.jl")
@everywhere LinearAlgebra.BLAS.set_num_threads(1)


function scan_gr(p::Params)
    @everywhere @unpack q, θ, E, R, N, seed = $p 
    rng = [MersenneTwister(seed + i) for i in 1:nprocs()]
    g_r = SharedArray{Float64}(length(E), N÷2)  
    g_r_sq = similar(g_r)
    @sync @distributed for i in 1:length(E)
            g_r[i, :], g_r_sq[i, :] = cor_tmm(θ = θ, E = E[i], N = N, q = q, R = R, rng = rng[myid()])
    end 
 
    return g_r, g_r_sq
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "config"
            help = "configuration file"
            required = true
        "--out"
            help = "Output file name"
            arg_type = String
            default = "gr.dat"
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()
    config = JSON.parsefile(args["config"])

    if haskey(config, "prec")
        setprecision(config["prec"])
    else 
        setprecision(250)
    end

    p = read_config(config)
    println(nprocs())
    @time g_r, _ = scan_gr(p)
    writedlm(args["out"], g_r)
end

main()
