using BunchFactorization
using LinearAlgebra
using Random

BLAS.set_num_threads(1)
#
println("________________________TEST________________________")
println()
rng = MersenneTwister(1234)
N = 3000
A = rand(rng, N, N)
A = A-A'
A_copy = copy(A)
A_copy2 = copy(A)

F  = bunch!(copy(A), pivot = true)
F2 = lu(A)

@time F  = bunch!(A_copy, pivot = true)
@time F2 = lu(A)
resid = A[F.p,F.p] - F.L*F.D*F.L'
resid2 = F2.P*A - F2.L*F2.U

b = rand(rng, N)
x = solve(F, b)
x2 = F2\b

println("N = 3000")
println("BUNCH Factorization Accuracy:  ", norm(resid)/norm(A))
println("LU    Factorization Accuracy:  ", norm(resid2)/norm(A))
println("Linear solver accuracy(BUNCH): ", norm(A*x - b)/norm(x))
println("Linear solver accuracy(LU):    ", norm(A*x2 - b)/norm(x))

N = 3001
A = rand(rng, N, N)
A = A-A'
A_copy = copy(A)
@time F  = bunch!(A_copy, pivot = true)

resid = A[F.p,F.p] - F.L*F.D*F.L'
println("norm: ", norm(resid)/norm(A))
