using BunchFactorization
using LinearAlgebra
using Random
using RandomSkewMatrices
using BandedMatrices
using MultiFloats

BLAS.set_num_threads(1)


#
println("________________________TEST________________________")
println()
rng = MersenneTwister(1234)
N = 10001
l = 6
F = Float64x3
A=F.(BandedMatrix(compact_chain(N, rng), (l^2, l^2)))
B = A[1:end-1, 1:end-1]
b = -A[1:end-1, end]
B_copy = copy(B)

F  = bunch!(copy(B),l, pivot = false)
x  = solve(F,l^2,b)
@time F  = bunch!(B_copy, l, pivot = false)
@time x = solve(F,l^2,  b)
push!(x, 1)

x = x/sqrt(sum(abs2, x))

resid2 = A*x

println("Linear solver accuracy: ", norm(resid2))
