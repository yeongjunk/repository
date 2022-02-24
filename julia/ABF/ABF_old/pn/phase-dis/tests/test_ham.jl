
using Plots
using Random
include("../library/abf2_ham.jl")
include("../library/abf2_pnscan.jl")
include("../library/abf2_cob.jl")

detangle = false
logscale = false
linscale_squared = true
θ = 0.25; N = 100; W = 0.01 ;i = 50;V = 0.1

r_arr = W*(rand(2N) .- 0.5)
H = create_ham_abf2(θ,θ,N)
add_diag!(H, r_arr)
H = convert(Array{ComplexF64,2}, H)
r_arr2 = V*(rand(2N, 3) .- 0.5)

add_hop_dis!(H,r_arr2)
eig = eigen(Hermitian(H));

a, b = splitpsi(eig.vectors[:,i])

if detangle == true
    U = unitary(θ)
    U = convert(Array{ComplexF64,2}, U)
    a, b = unitary(a, b,U')
    a, b = uc_redef_n(a, b)
    a, b = unitary(a, b, U')
end

if logscale == true
    plt = scatter(log10.(abs.(a)), label = "a", color = :blue)
    scatter!(plt,log10.(abs.(b)), label = "b", color = :red)
else
    if linscale_squared == true
        plt = scatter(abs.(a).^2, label = "a", color = :blue)
        scatter!(plt,abs.(b).^2, label = "b", color = :red)
    else
        plt = scatter(a, label = "a", color = :blue)
        scatter!(plt,b, label = "b", color = :red)
    end
end
title!("E = $(round(eig.values[i], digits = 6))")
display(plt)
plt2 = plot()
# plot!(plt2, eig.values)
# display(plt2)
testU = isapprox(U*U', I(2))
testH = ishermitian(H)
println("Is U times U_adj I?: $(testU)")
println("Is H hermitian?: $(testH)")
println("E = $(eig.values[i])")
println("PN = $(compute_pn((eig.vectors[:,i])))")

println(minimum(eig.values))
