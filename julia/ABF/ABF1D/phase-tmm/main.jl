using Distributed
using JSON, ArgParse
using CSV
using DataFrames
using SharedArrays
@everywhere include("./xi.jl")
@everywhere include("./xi_estimate.jl")
@everywhere LinearAlgebra.BLAS.set_num_threads(1)


function scan_xi(p::Params)
    @everywhere @unpack θ, E, R, N, seed = $p 
    rng = [MersenneTwister(seed + i) for i in 1:nprocs()]
    xi = SharedArray{Float64}(length(E), R)  
    for i in 1:length(E)
       @sync @distributed for r in 1:R
            xi[i, r] = compute_xi(θ = θ, E = E[i], N = N[i], rng = rng[myid()])
       end 
    end 
    return mean(xi, dims = 2), std(xi, dims=2)
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
            default = "xi.csv"
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()
    config = JSON.parsefile(args["config"])
    p = read_config(config)
    if p.N_auto == true
        m, a = estimate_xi_coefs(θ = p.θ, N = 10^6)
        p.N .= estimate_N.(p.E, a, m)
        idx = findall(x -> x < 10^6, p.N)
        p.N[idx] .= 10^6
        println(p.N)
    end
    println(nprocs())
    @time xi, xi_std = scan_xi(p)
    CSV.write(args["out"], DataFrame(E = p.E, xi = vec(xi), xi_std = vec(xi_std), xi_ste = vec(xi_std)/sqrt(p.R)))
end

main()
