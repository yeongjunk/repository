using LinearAlgebra
using StatsBase

include("./xi.jl")
LinearAlgebra.BLAS.set_num_threads(1)

#Compute rough estimate of xi = a*E^b
function estimate_xi_coefs(; θ = 0.25, E1 = 10^-4, E2 = 10^-5, N = 10^7, R = 2)
    xi_1 = mean([compute_xi(θ = θ, E = E1, N = N) for i in 1:R])
    xi_2 = mean([compute_xi(θ = θ, E = E2, N = N) for i in 1:R])
    log_xi = log10.([xi_1; xi_2])
    log_E = [-4.; -5.] 
    Y = log_xi
    X = zeros(2, 2); X[:, 1] = log_E;X[:, 2] .= 1. 
    m, a = X\Y # m: power, a = prefactor
    a = 10^a
return m, a
end

estimate_xi(E, a, b) = a*E^b
estimate_N(E, a, b) = trunc(Int64, 10^6*estimate_xi(E, a, b)) # Roughly 10^6 times the localization length is a good N to get error within 10%?

