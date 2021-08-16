using JLD
using JSON
using ArgParse
include("library/abf2_pnscan.jl")

function main(args)
    opts = ArgParseSettings(description="Scan and compute pn for all parameters of nu=2 ABF")
    @add_arg_table! opts begin
    "c"
        help = "configuration"
        arg_type = AbstractString
    end

    # Parse the arguments
    popts   = parse_args(opts)
    config  = JSON.parsefile(popts["c"])

    BLAS.set_num_threads(1) #Turn off BLAS multi-threading

    full_p = readconfig(config)   # create the Parameters struct from Dict config

    V, θ = expand_params(full_p)
    @time for i in 1:full_p.V_num, j in 1:full_p.θ_num
            p = full_to_one(full_p, i,j)
            fn_out = "V_th_$(i)_$(j).jld"
            df = scan_pn(p)
            dict = Dict("W" => p.W, "V" => V[i], "θ" => θ[j], "data" => df)
            JLD.save(fn_out, dict)
    end

    println("DATA SAVED")
end

main(ARGS)
