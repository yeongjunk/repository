using LinearAlgebra
using StatsBase

include("./xi.jl")
LinearAlgebra.BLAS.set_num_threads(1)

#Compute rough estimate of xi = a*E^b
function estimate_xi(; θ = 0.25, E = [0.01], N = 10^5, q = 4)
   xi = similar(E)
   for i in 1:length(E)
       xi[i] = compute_xi(θ = Float64(θ), E = E[i], N = N,  q = q)
   end
   return xi
end

function estimate_N(; θ = 0.25, E = 0.01, N = 10^5, q = 4) 
   Ns = trunc.(Int, 10^5*estimate_xi(θ = θ, E = E, N = N, q = q))
   idx = findall(x -> x < 1000000, Ns)
   Ns[idx] .= 1000000
   return Ns
end
