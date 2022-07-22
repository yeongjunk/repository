# Diagonalize scale free model of phase disordered d = 1 nu = 2 ABF using Lanczos algorithm.
# and compute GIPR(generalized IPR) of moment q and box size l.

using LinearAlgebra, SparseArrays, Arpack
using LinearMaps
using Random
using DataFrames, CSV
using ArgParse, JSON
# Custom modules
using ABF
using Lattice
using PN
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())

ENV["JULIA_COPY_STACKS"] = 1
include("./ed.jl") # read parameters from configuration file

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
    abf_scan(p)
end

main(ARGS)
