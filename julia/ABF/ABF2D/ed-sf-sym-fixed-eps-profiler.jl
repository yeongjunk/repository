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

include("./ed-sf-sym-fixed-eps-func.jl") # read parameters from configuration file
config  = JSON.parsefile("ed-sf-sym-fixed-eps-config")
p = readconfig(config)
@profiler abf3d_scan(p)
