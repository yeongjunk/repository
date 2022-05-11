using JSON
using MAT
using ArgParse

include("./xi.jl")
include("./params.jl")

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
    E = range(config["E"][1], config["E"][2], length = config["E"][3])
    E = 10 .^E
    R = config["R"]
    xi = Array{Float64}(undef, length(E), length(R))
    fill!(xi, 0.)
    for r in 1:config["R"]    
        @Threads.threads for i in 1:length(E)
            p = Params(E[i], Float64(config["th"]), Float64(config["V"]), Float64(config["W"]), Int64(config["q"]), Int64(config["N"]), Int64(config["seed"]+r))
            xi[i,r] = loc_length(p)
        end
    end
    matwrite("./xi.mat", Dict("config" => config, "xi" => xi, "E" => collect(E)))
    println("Scan finished. Saved in $(pwd())")
end

main(ARGS)

