using PN
using Plots
include("./ed.jl")
include("./rgf.jl")

## CHAIN
N = 500; W = 1.
H = ham_chain(N = N)
D = spdiagm(0 => W*(rand(size(H,1)) .- 0.5))
E, ψ = eigen(Hermitian(Matrix(H + D)))
pn = compute_pns(ψ)

ξ = similar(E)
for i in 1:length(E)
    ξ[i] = rgf_strip.(E = E[i], M = 1, N = 1E5)
end

p_chain = scatter(E, pn, label =  "ED")
scatter!(p_chain, E, ξ, label = "RGF")
## STRIP
#------------------------ ED ------------------------#
M = 3; N = 400; W = 2.
H = ham_strip(M = M, N = N)
D = spdiagm(0 => W*(rand(size(H,1)) .- 0.5))
E, ψ = eigen(Hermitian(Matrix(H + D)))
pn = compute_pns(ψ)
#------------------------ RGF ------------------------#
E_rgf = -4:0.1:4
ξ = similar(E_rgf)
for i in 1:length(E_rgf)
    ξ[i] = rgf_strip(M = M, W = 1., E = E_rgf[i], N = 1E5)
end

p_strip = scatter(E, pn)
plot!(p_strip, E_rgf, ξ, lw = 2)
plot!(p_strip, E_rgf, M*ξ, lw = 2)

## Bar
M = 3; N = 400; W = 5.
H = ham_bar(M = M, N = N)
D = spdiagm(0 => W*(rand(size(H,1)) .- 0.5))
E, ψ = eigen(Hermitian(Matrix(H + D)))
pn = compute_pns(ψ)
#------------------------ RGF ------------------------#
E_rgf = -4:0.1:6
ξ = similar(E_rgf)
for i in 1:length(E_rgf)
    ξ[i] = rgf_bar(M = M, W = 5., E = E_rgf[i], N = 1E5)
end

p_bar = scatter(E, pn, c = :white)
plot!(p_bar, E_rgf, ξ, lw = 3, c = :red)
plot!(p_bar, E_rgf, M^2*ξ, lw = 3, c = :blue)
