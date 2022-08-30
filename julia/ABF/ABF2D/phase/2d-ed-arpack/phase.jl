using Lattices
using ABF
using Random
using SparseArrays
using LinearAlgebra  
using Distributions

function phase_dis(H; V = 1., rng = nothing)
    F = eltype(H)
    D = convert.(Complex{F}, H)
    if rng == nothing
        rng = MersenneTwister()
    end 
    rows = rowvals(D)
    vals = nonzeros(D)
    m, n = size(D)
    for j = 1:n 
       for i in nzrange(D, j)
          row = rows[i]
          # println("$row ,", "$j")
          if row > j 
              vals[i] = im*V*vals[i]*(rand(rng, F) .- 0.5)
          elseif row <= j
              vals[i] = F(0.)
          end 
       end 
    end 
    return D + D'
end


function ham_sf(;L::Integer = 10, V = 1, θ::F = 0.25,rng = Random.GLOBAL_RNG) where F <: AbstractFloat
    ltc = Lattice2D(L, L, 2)
    H, U = ham_fe(ltc, F(-2.), F(0.), θ)
    D = phase_dis(H, V = V, rng = rng)
    H_sf = project(U'*D*U)
    return H_sf
end
