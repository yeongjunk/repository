using Bunch
using LinearAlgebra

N = 1000
A = rand(N, N)
A = A-A'
A_copy = copy(A)
p, L, D = bunch!(A_copy, pivot=true)

N = 3000
A = rand(N, N)
A = A-A'
A_copy = copy(A)
@time p, L, D = bunch!(A_copy, pivot = true)

resid = A[p,p] - L*D*L'
println("norm: ", norm(resid)/norm(A))
