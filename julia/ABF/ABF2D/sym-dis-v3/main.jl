using LinearAlgebra, SparseArrays 
using LinearMaps
using Random
using DataFrames, CSV
using ArgParse, JSON
using JLD
# Custom modules
using ABFSym
using Lattice
using PN
ENV["JULIA_COPY_STACKS"] = 1

include("./lib.jl") # read parameters from configuration file
include("./gen-e-params.jl")

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

    LinearAlgebra.BLAS.set_num_threads(config["num_blas"])

    fn_e_params = "./e-params.jld" 
    if config["E_param_gen"] == true
        println("Generating energy parameters.")
        @time E_bw, E_c, E_del = gen_e_params(p)
        JLD.save(fn_e_params, "E_bw" ,E_bw, "E_c",  E_c, "E_del" ,E_del)
        display(E_c)
    else
        if config["auto_E"] == true
            println("Reading energy parameters")
            #---- Read automatically generated energy parameters ----#
            eparams = JLD.load(fn_e_params)  
            E_bw =  eparams["E_bw"]
            E_c =   eparams["E_c"]
            E_del = eparams["E_del"]

            p_E = EParams(E_bw, E_c, E_del)

            @time abf2d_scan(p, p_E)
        else 
            @time abf2d_scan(p)
        end
    end
end

main(ARGS)
