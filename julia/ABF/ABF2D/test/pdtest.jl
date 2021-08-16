using LinearAlgebra
using Plots
using SparseArrays
using Random
using DataFrames
using StatsBase
using Arpack

gr(fmt = :png)
L = 64
l = Lattice2D(L, L, 2)
l_p = Lattice2D(L, L, 1)
W = 0.
V = 1.

@profiler H, U = ham_fe(l, -1, 1, 0.25)
H = convert.(ComplexF64, H)
add_pd!(H, V, rng)
H .= U'*H*U

H_p = project(H)
droptol!(H_p, 1E-12)
