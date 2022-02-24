using LinearAlgebra
using Plots


ϕ = 0.1; N = 100
include("../library/abf2_ham.jl")
U = unitary_block(pi/4)
T = uc_redef(N)
U_full = blockdiagonal(U, N)
t = exp(2π*im*ϕ)

H_det_block = Array{ComplexF64}([0 -t -conj(t); -conj(t) 0 -t; -t -conj(t) 0])
H = blockdiagonal(H_det_block, N)
H = T*H*T'
H = U_full*H*U_full'

# H = H .+ Diagonal(W*(rand(300) .- 0.5))
val, vect = eigen(Hermitian(Array(H)))

p = plot(val)
