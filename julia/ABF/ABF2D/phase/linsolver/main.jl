using Distributed
using JSON, ArgParse
using CSV, DataFrames
using SharedArrays
using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)
@everywhere using LinearAlgebra
@everywhere using MultiFloats
@everywhere MultiFloats.use_bigfloat_transcendentals()
x = LOAD_PATH
@everywhere lp = @fetchfrom 1 x 
@everywhere append!(LOAD_PATH, lp)
@everywhere include("./gr.jl")
@everywhere include("./phase.jl")

function grscan(;L = 20, θ = 0.25, R = 10, rng=Random.GLOBAL_RNG)
    grs = Array{Float64}[]
    ltc = Lattice2D(L, L, 1)
    @time for r in 1:R
        psi = compute_psi(L = L, θ = θ, rng = rng);
        lw = LatticeWave(ltc, log.(abs.(psi)));
        g_r = filter(!isnan, eig_corr_full(lw));
        push!(grs, g_r)
    end
    return vec(mean(reduce(hcat, grs), dims = 2))
end


@everywhere function grscan_sum(;L = 20, θ = 0.25, R = 10, rng=Random.GLOBAL_RNG)
    grs = Array{Float64}[]
    ltc = Lattice2D(L, L, 1)
    @time for r in 1:R
        psi = compute_psi(L = L, θ = θ, rng = rng);
        lw = LatticeWave(ltc, log.(abs.(psi)));
        g_r = filter(!isnan, eig_corr_full(lw));
        push!(grs, g_r)
    end
    return vec(sum(reduce(hcat, grs), dims = 2))
end

function grscan_distributed(p::Params)
    @everywhere @unpack θ, R, L, seed = $p 
    rng = [MersenneTwister(seed + i) for i in 1:nprocs()]
    R_dist = R ÷ nprocs()
    futures = Future[]
    for i in nprocs():-1:1
        push!(futures, remotecall(grscan_sum ,i ,L=L, θ=θ, R=R_dist, rng = rng[i]))
    end  
    grs = Array{Float64}[]
    for i in 1:nprocs()
        push!(grs, fetch(futures[i]))
    end  
    grs = reduce(hcat, grs)
    g_r = vec(sum(grs, dims = 2))/(R_dist*nprocs()) 
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
            default = "gr.csv"
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
    println("Number of processors: ",nprocs())
    @time g_r = grscan_distributed(p)
    CSV.write(args["out"], DataFrame(r = vec(1:length(g_r)), g_r =  g_r))
end

main()
