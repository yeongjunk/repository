const modulepath =  "/Users/pcs/codes/chain/ABF2D/module"

using Plots
using LinearAlgebra
using SparseArrays
using Random
using DataFrames
plotly()

push!(LOAD_PATH, modulepath)
using ABF2D
using Lattice
using PN

function abf_pn!(l_fe::Lattice2D, l_prj::Lattice2D, U,V, rng; q = 2)
    T = phase_dis(l_fe, V, rng)
    H_fe = copy(l_fe.H)
    H_fe .+= T
    H_fe .= U'*H_fe*U   #detangle
    project!(l_prj.H, H_fe)
    E, PN = eig_pn!(Hermitian(Array(l_prj.H)), q = q)
    df = DataFrame()
    df.E = E; df.PN = PN
    return df
end


L = 41; θ = 0.25; rng = MersenneTwister(1234)

l = Lattice2D(L, L, 2)
ham_fd!(l, -1, 1)
U = U_fe(l, 0.25)
l.H .= U*l.H*U'

l_prj = Lattice2D(L,L,1)

df = abf_pn!(l, l_prj,U, 1., rng, q = 2.)
@time df = abf_pn!(l, l_prj,U, 0.001, rng)
@profiler df = abf_pn!(l, l_prj,U, 0.001, rng)

scatter(df.E, df.PN, markersize = 1, xlabel = "E", ylabel = "PN")
