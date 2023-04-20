using BunchFactorization
using LinearAlgebra
using Random
using RandomSkewMatrices
BLAS.set_num_threads(1)
#
println("________________________TEST________________________")
println()
rng = MersenneTwister()
N = 1001

A = Matrix(compact_chain(N, rng))

# A = rand(N, N) 
# A = A - A' 
B = A[1:end-1, 1:end-1]
b = -A[1:end-1, end]
A_copy = copy(B)
A_copy2 = copy(B)

F  = bunch!(copy(A), pivot = true)
F2 = lu(A)
F3  = bunch2!(copy(A), pivot = true)

@time F  = bunch!(A_copy, pivot = true)
@time F2 = lu(B)
@time F3 = bunch2!(A_copy2, pivot= true)
resid = B[F.p,F.p] - F.L*F.D*F.L'
resid2 = F2.P*B - F2.L*F2.U
resid3 = B[F3.p,F3.p] - F3.L*F3.D*F3.L'

x = solve(F, b)
x2 = F2\b
x3 =solve(F3, b) 
push!(x, 1)
push!(x2, 1)
push!(x3, 1)

println("N = 3000")
println("BUNCH Factorization Accuracy:  ", norm(resid)/norm(B))
println("LU    Factorization Accuracy:  ", norm(resid2)/norm(B))
println("BUNCH+ROOK    Factorization Accuracy:  ", norm(resid3)/norm(B))

println("Linear solver accuracy(BUNCH): ", norm(A*x)/norm(x))
println("Linear solver accuracy(LU):    ", norm(A*x2)/norm(x2))
println("Linear solver accuracy(BUNCH+ROOK):    ", norm(A*x3)/norm(x3))

println("N = 3001")
N = 3001
A = rand(rng, N, N)
A = A-A'
A_copy = copy(A)
@time F  = bunch!(A_copy, pivot = true)
resid = A[F.p,F.p] - F.L*F.D*F.L'
println("norm: ", norm(resid)/norm(A))

