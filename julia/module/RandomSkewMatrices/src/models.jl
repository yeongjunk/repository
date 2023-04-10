using Lattices 
using SparseArrays
using Distributions
using Random

export compact_chain, erbm, prbm, quasiperiodic_compact_chain

function compact_chain(N, rng;  dist::Distribution=Uniform(-0.5, 0.5), l=2, bc=:open)
    ltc = Lattice1D(N, 1)
    I = Int64[]; J=Int64[]; V=Float64[]
    for n in 1:ltc.N    
       for i in 1:l 
           n+i > ltc.N && bc == :open && break
           push!(I, index(ltc, (n, 1))) 
           push!(J, index(ltc, (n+i, 1))) 
           push!(V, rand(rng, dist))
       end
    end
    H = sparse(I, J, V, ltc.N, ltc.N)
    return H - H' 
end


function erbm(N::Integer, rng; dist::Distribution=Uniform(-0.5, 0.5), exponent = 1.0)
    ltc = Lattice1D(N, 1)
    I = Int64[]; J=Int64[]; V=Float64[]
    for n in 1:ltc.N    
       for i in 1:ltc.N 
           n+i > ltc.N && continue
           push!(I, index(ltc, (n, 1))) 
           push!(J, index(ltc, (n+i, 1))) 
           push!(V, exp(-exponent*(i-1))*rand(rng, dist))
       end
    end
    H = sparse(I, J, V, ltc.N, ltc.N)
    return H - H' 
end

function prbm(N::Integer, rng::AbstractRNG; dist::Distribution=Uniform(-0.5, 0.5), exponent = 2.0) 
    ltc = Lattice1D(N, 1)
    I = Int64[]; J=Int64[]; V=Float64[]
    for n in 1:ltc.N    
       for i in 1:ltc.N 
           n+i > ltc.N && continue
           push!(I, index(ltc, (n, 1))) 
           push!(J, index(ltc, (n+i, 1))) 
           push!(V, (i)^(-exponent)*rand(rng, dist))
       end
    end
    H = sparse(I, J, V, ltc.N, ltc.N)
    return H - H' 
end

function quasiperiodic_compact_chain(N::Integer, rng::AbstractRNG; a = (sqrt(5) - 1)/2, phi = 0.0)
    ltc = Lattice1D(N, 1)
    I = Int64[]; J=Int64[]; V=Float64[]
    randphase = 2pi*rand(rng)
    for n in 1:ltc.N    
       for i in 1:2 
           n+i > ltc.N && break
           push!(I, index(ltc, (n, 1))) 
           push!(J, index(ltc, (n+i, 1))) 
           push!(V, cos(a*pi*n+i*phi+randphase))
       end
    end
    H = sparse(I, J, V, ltc.N, ltc.N)
    return H - H' 
end
