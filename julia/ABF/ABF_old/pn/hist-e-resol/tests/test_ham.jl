
using Plots
using Random
include("../library/abf2_ham.jl")
include("../library/abf2_pnscan.jl")
include("../library/abf2_cob.jl")

seed = rand(1:100000);
detangle = true
logscale = true
linscale_squared = true
θ = 0.25; N = 100; W = 10^-10;i = 25
# almost CLS
seed = 97560; i = 50; θ = 0.25; N = 100; W = 10^-10
#seed = 32229; i = 50; θ = 0.25; N = 100; W = 10^-10
#seed = 70652; i = 50; θ = 0.25; N = 100; W = 10^-10


W = convert(Float64, W)
H = create_ham_abf2(θ,θ,N)

rng = MersenneTwister(seed)
r_arr = W*(rand(rng, 2N) .- 0.5)

H_dis = add_diag(H,r_arr)
eig = eigen(Symmetric(H_dis));

a, b = splitpsi(eig.vectors[:,i])

if detangle == true
    U = unitary(θ)
    a, b = unitary(a, b,U')
    a, b = uc_redef_n(a, b)
    a, b = unitary(a, b, U')
end

if logscale == true
    plt = scatter(log10.(abs.(a)), label = "a", color = :blue)
    scatter!(plt,log10.(abs.(b)), label = "b", color = :red)
else
    if linscale_squared == true
        plt = scatter(a.^2, label = "a", color = :blue)
        scatter!(plt,b.^2, label = "b", color = :red)
    else
        plt = scatter(a, label = "a", color = :blue)
        scatter!(plt,b, label = "b", color = :red)
    end
end
title!("E = $(round(eig.values[i], digits = 6))")
display(plt)

testU = isapprox(U*U', I(2))
testH = ishermitian(H)
println("Is U times U_adj I?: $(testU)")
println("Is H hermitian?: $(testH)")
println("E = $(eig.values[i])")
println("PN = $(compute_pn(eig.vectors[:,i]))")
println("seed = $(seed)")
