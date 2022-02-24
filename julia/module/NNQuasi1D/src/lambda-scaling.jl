using Plots
using LinearAlgebra
include("./ed.jl")
include("./rgf.jl")

M = [2 4 8 16 32]
W = vec(1.:0.2:3.)

ξ = Array{Float64}(undef, length(W), length(M))

for i in 1:length(M), j in 1:length(W)
    ξ[j, i] = rgf_strip(M = M[i], W = W[j], N = 1E5)
end

p = scatter()
for i in 1:length(W)
    scatter!(p, vec(M), ξ[i,:]./vec(M), yaxis =:log10, xaxis = :log10)
end

p
