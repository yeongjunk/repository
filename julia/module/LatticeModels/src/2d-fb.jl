using Lattices
using SparseArrays

export kagome 

"""
keyword arguments
t1: main hoppings
t2: auxiliary hopping that controls the FB energy
e_a: onsite energy at a(middle sublttice)
e_bc: onsite energy at b, c (upper, lower sublattice)
"""
function kagome(;N::Integer=10, M::Integer=10, t::F = 1.0) where F
    ltc = Lattice2D(N, M, 3)
    I = Int64[]; J=Int64[]; V=F[]
    a, b, c = 1, 2, 3

    for n in 1:ltc.N, m in 1:ltc.M    
       n1 = mod1(n+1, ltc.N)
       m1 = mod1(m+1, ltc.M)

       push!(I, index(ltc, (n,m,a))) 
       push!(J, index(ltc, (n,m,b))) 
       push!(V, -t) 

       push!(I, index(ltc, (n,m,a))) 
       push!(J, index(ltc, (n,m,c))) 
       push!(V, -t)

       push!(I, index(ltc, (n,m,b))) 
       push!(J, index(ltc, (n,m,c))) 
       push!(V, -t)


       push!(I, index(ltc, (n,m,c))) 
       push!(J, index(ltc, (n1,m,a))) 
       push!(V, -t) 

       push!(I, index(ltc, (n,m,b))) 
       push!(J, index(ltc, (n,m1,c))) 
       push!(V, -t) 

       push!(I, index(ltc, (n,m,b))) 
       push!(J, index(ltc, (n1,m1,a))) 
       push!(V, -t) 
    end
    H = sparse(I, J, V, ltc.N*ltc.M*ltc.U, ltc.N*ltc.M*ltc.U)
    return H + H'
end
