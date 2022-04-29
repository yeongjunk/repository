using Distributed
using JSON, ArgParse
using MAT 
using DataFrames
using SharedArrays
@everywhere include("./xi.jl")
@everywhere LinearAlgebra.BLAS.set_num_threads(1)


function scan_lyap(p::Params)
    @everywhere @unpack q, θ, E, R, N, seed = $p 
    rng = [MersenneTwister(seed + i) for i in 1:nprocs()]
    lyap = SharedArray{Float64}(length(E), length(N), R)  
    for i in 1:length(E)
        for j in 1:length(N) 
           @sync @distributed for r in 1:R
               lyap[i, j, r] = compute_lyap(θ = θ, E = E[i], N = N[j], q = q, rng = rng[myid()])
            end 
        end
    end 
    return lyap 
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
            default = "lyap.mat"
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()
    config = JSON.parsefile(args["config"])
    if haskey(config, "prec")
        setprecision(config["prec"])
    end

    p = read_config(config)
    println(nprocs())
    @time lyap = scan_lyap(p)
    lyap = Array(lyap)
    dct = Dict("E" => Float64.(p.E), "N" => p.N, "R" => p.R, "q" => p.q, "lyap" => lyap)
    display(dct)
    matwrite(args["out"], dct)  
end

main()
