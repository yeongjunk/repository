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
N = rand(20:3001) 
# N = 17 
l = rand(1:9) 
# l = 1 
println("N: ", N, ", l: ", l)
F = Float64
A=F.(BandedMatrix(compact_chain(N, rng, l = l), (l+1, l+1)))
A_copy = copy(A)

bunch!(copy(A), l, pivot = false)
@time F  = bunch!(A_copy, l, pivot = false)
println("accuarcy(nopivot): ", norm(A - F.L*F.D*F.L'))



A= Matrix(compact_chain(N, rng, l = l)) 
A_copy = copy(A)

F  = bunch!(copy(A), l, pivot = :partial)
@time F  = bunch!(A_copy, l, pivot = :partial)
println("accuarcy(partial): ", norm(A[F.p, F.p] - F.L*F.D*F.L'))
# display(F.L)

A= Matrix(compact_chain(N, rng, l = l)) 
A_copy = copy(A)

F  = bunch!(copy(A), l, pivot = :partial2)
@time F  = bunch!(A_copy, l, pivot = :partial2)
println("accuarcy(partial): ", norm(A[F.p, F.p] - F.L*F.D*F.L'))
