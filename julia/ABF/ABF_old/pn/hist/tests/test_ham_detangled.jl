using Plots
gr()
using Random
using LinearAlgebra
include("../library/abf2_ham.jl")
include("../library/abf2_cob.jl")
θ = 1/4;
N = 100;
W = 0;
W = convert(Float64, W)
H0 = [1.0 0.0; 0.0 -1]

H = create_ham_abf2(θ, θ, N)
display("text/plain", H)
U = unitary(θ)

testU = isapprox(U * U', I(2))
testH = ishermitian(H)

println("Is U times U_adj I?:" * " $(testU)")
println("Is H hermitian?:" * " $(testH)")


rng = MersenneTwister(1234)

r_arr = W * (rand(rng, 2N) .- 0.5)
H_dis = add_diag(H, r_arr)

eig = eigen(Symmetric(H_dis));

psi = eig.vectors[:, 3]
U_adj = unitary(-θ)
a, b = splitpsi(psi)
p = scatter(a, legend = false)
scatter!(p, b, legend = false)

a, b = unitary(a, b, U_adj)
a, b = uc_redef_n(a, b)
a, b = unitary(a, b, U_adj)

p1 = scatter(a, legend = false)
scatter!(p1, b, legend = false)
