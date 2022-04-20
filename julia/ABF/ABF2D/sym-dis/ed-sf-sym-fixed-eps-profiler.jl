using LinearAlgebra, SparseArrays, Arpack
using LinearMaps
using Random
using DataFrames, CSV
using ArgParse, JSON
# Custom modules
using ABFSym
using Lattice
using PN
LinearAlgebra.BLAS.set_num_threads(1)
rdir = @__DIR__
include("./ed-sf-sym-fixed-eps-func-threaded-search.jl") # read parameters from configuration file
config  = JSON.parsefile(rdir*"/sym-config2")
p = readconfig(config)

@profiler abf3d_scan(p)
