using LinearAlgebra, SparseArrays, Arpack
using LinearMaps
using Random
using DataFrames, CSV
using ArgParse, JSON
# Custom modules
using ABFSym
using Lattice
using PN
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
ENV["JULIA_COPY_STACKS"] = 1
include("./ed-sf-sym-fixed-eps-func.jl") # read parameters from configuration file
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
    abf3d_scan(p)
end

main(ARGS)
