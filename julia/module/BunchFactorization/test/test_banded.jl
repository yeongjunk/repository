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
rng = MersenneTwister()
N = 35
l = 4
F = Float64x3
A=F.(BandedMatrix(compact_chain(N, rng, l = l), (l+1, l+1)))
B = A[1:end-1, 1:end-1]
b = -A[1:end-1, end]
B_copy = copy(B)
F  = bunch!(copy(B),l, pivot = false)
display(F.L)

x  = solve(F,l+1, b)

@time F  = bunch!(B_copy, l, pivot = false)
@time x = solve(F,l+1,b)
push!(x, 1)

x = x/sqrt(sum(abs2, x))

resid2 = A*x

println("Linear solver accuracy: ", norm(resid2))
