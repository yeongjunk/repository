using Distributed
using JSON, ArgParse
using CSV
using DataFrames
using SharedArrays
using DelimitedFiles
@everywhere using LinearAlgebra
@everywhere using MultiFloats
@everywhere MultiFloats.use_bigfloat_transcendentals()
x = LOAD_PATH
@everywhere lp = @fetchfrom 1 x 
@everywhere append!(LOAD_PATH, lp)
@everywhere include("./gr.jl")

function scan_gr(p::Params)
    @everywhere @unpack θ, R, L, seed = $p 
    rng = [MersenneTwister(seed + i) for i in 1:nprocs()]
    g_r = SharedArray{Float64}(L÷2)  
    g_r_sq = similar(g_r)
    if length(rng) == 1
        g_r, _ = cor_tmm(θ = θ, L = L, R = R, rng = rng[1])
    else
        g_r = cor_tmm(θ = θ, L = L, R = R, rng = rng)
    end 

    return g_r
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
    @time g_r = scan_gr(p)
    writedlm(args["out"], g_r)
end

main()
