using Lattice
using Random
using SparseArrays
using LinearAlgebra  
using Distributions

function nnhopping(p::AbstractArray{T}, s::T, c::T, s4::T,c4::T, s2c2::T) where T <: Real 
    t_1 = 2im*(p[1, 1]*s4 - p[1, 1]*c4 + p[2, 1]*s2c2 + p[2, 2]*s2c2 
        + p[3, 1]*c4 - p[3, 2]*s2c2 - p[4, 1]*s2c2 + p[4, 2]*s4
        - p[5, 1]*c4 - p[5, 2]*s4)*s2c2
    return t_1
end

function nnnhopping(p::AbstractArray{T}, s::T, c::T, s4::T,c4::T, s2c2::T) where T <: Real 
    t_2 = 2im*(-p[2, 1] + p[3, 1] + p[4, 1] - p[5, 1])*s4*c4
    return t_2 
end


function ham_sf_obc(;L::Integer = 10, θ::F = 0.25,rng = Random.GLOBAL_RNG) where F <: AbstractFloat
    dist = Uniform(F(-0.5), F(0.5)) 
    # pre-calculation of sines and cosines for performance optimization 
    s = sinpi(θ); c = cospi(θ) 
    s4 = s^4; c4 = c^4; s2c2 = s^2*c^2 
    # Initial hoppings 
    p = Array{F}(undef, 5, 2) 
    rand!(rng, dist, p)

    ltc = Lattice1D(L, 1)
    I, J = Int64[], Int64[]
    V = Complex{F}[]
    for l in 1:L-2
        push!(I, index(ltc, (l, 1)))
        push!(J, index(ltc, (l+1, 1)))
        push!(V, nnhopping(p, s, c, s4, c4, s2c2))   
        
        push!(I, index(ltc, (l, 1))) 
        push!(J, index(ltc, (l+2,1)))
        push!(V, nnnhopping(p, s, c, s4, c4, s2c2))   
        p[:, 2] = @view p[:, 1]
        rand!(rng, dist, @view p[:, 1])
        if l == L
            p[:, 1] = p_final    
        end
    end
    l = L-1
        push!(I, index(ltc, (l, 1)))
        push!(J, index(ltc, (l+1, 1)))
        push!(V, nnhopping(p, s, c, s4, c4, s2c2))   

    H = sparse(I, J, V, L, L) 
    return (H + H')/2
end

function ham_sf_pbc(;L::Integer = 10, θ::F = 0.25,rng = Random.GLOBAL_RNG) where F <: AbstractFloat
    # Create uniform distribution pdf
    dist = Uniform(F(-0.5), F(0.5)) 
    # pre-calculation of sines and cosines for performance optimization 
    s = sinpi(θ); c = cospi(θ) 
    s4 = s^4; c4 = c^4; s2c2 = s^2*c^2 
    #hoppings 
    p = Array{F}(undef, 5, 2) 
    rand!(rng, dist, p)
    p_final = p[:, 1]

    ltc = Lattice1D(L, 1)
    I, J = Int64[], Int64[]
    V = Complex{F}[]
    for l in 1:L
        push!(I, index(ltc, (l, 1)))
        push!(J, index(ltc, (l+1, 1)))
        push!(V, nnhopping(p, s, c, s4, c4, s2c2))   
        
        push!(I, index(ltc, (l, 1))) 
        push!(J, index(ltc, (l+2,1)))
        push!(V, nnnhopping(p, s, c, s4, c4, s2c2))   
        p[:, 2] = @view p[:, 1] # hoppings are correlated so they have to be propagated to the next loop
        rand!(rng, dist, @view p[:, 1]) # new hoppings
        if l == L
            p[:, 1] = p_final    
        end
    end
    H = sparse(I, J, V, L, L) 
    return (H + H')/2
end
