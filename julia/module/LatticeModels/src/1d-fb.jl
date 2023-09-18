using Lattices
using SparseArrays

export diamond, stub

"""
keyword arguments
t1: main hoppings
t2: auxiliary hopping that controls the FB energy
e_a: onsite energy at a(middle sublttice)
e_bc: onsite energy at b, c (upper, lower sublattice)
"""
function diamond(;N::Integer=10, t1::F=1.0, t2::F=1.0, e_a::F = 1.0, e_bc::F = 1.0) where F
    ltc = Lattice1D(N, 3)
    I = Int64[]; J=Int64[]; V=F[]
    a, b, c = 1, 2, 3

    for n in 1:ltc.N    
       n1 = mod1(n+1, ltc.N)

       push!(I, index(ltc, (n, a))) 
       push!(J, index(ltc, (n, a))) 
       push!(V, e_a/2) 

       push!(I, index(ltc, (n, b))) 
       push!(J, index(ltc, (n, b))) 
       push!(V, e_bc/2) 
       push!(I, index(ltc, (n, c))) 
       push!(J, index(ltc, (n, c))) 
       push!(V, e_bc/2) 

       push!(I, index(ltc, (n, a))) 
       push!(J, index(ltc, (n, b))) 
       push!(V, t1) 

       push!(I, index(ltc, (n, a))) 
       push!(J, index(ltc, (n, c))) 
       push!(V, t1)

       push!(I, index(ltc, (n, b))) 
       push!(J, index(ltc, (n, c))) 
       push!(V, t2)

       push!(I, index(ltc, (n, b))) 
       push!(J, index(ltc, (n1, a))) 
       push!(V, t1) 

       push!(I, index(ltc, (n, c))) 
       push!(J, index(ltc, (n1, a))) 
       push!(V, t1) 
    end
    H = sparse(I, J, V, ltc.N*ltc.U, ltc.N*ltc.U)
    return H + H'
end


"""
keyword arguments
t1: main hoppings
t2: auxiliary hopping that controls the FB energy
e_a: onsite energy at a(middle sublttice)
e_bc: onsite energy at b, c (upper, lower sublattice)
"""
function stub(;N::Integer=10, t1::F=1.0, t2::F=1.0, e_ac::F = 1.0, e_b::F = 1.0) where F
    ltc = Lattice1D(N, 3)
    I = Int64[]; J=Int64[]; V=F[]
    a, b, c = 1, 2, 3

    for n in 1:ltc.N    
       n1 = mod1(n+1, ltc.N)

       push!(I, index(ltc, (n, a))) 
       push!(J, index(ltc, (n, a))) 
       push!(V, e_ac/2) 

       push!(I, index(ltc, (n, b))) 
       push!(J, index(ltc, (n, b))) 
       push!(V, e_b/2) 

       push!(I, index(ltc, (n, c))) 
       push!(J, index(ltc, (n, c))) 
       push!(V, e_ac/2) 

       push!(I, index(ltc, (n, a))) 
       push!(J, index(ltc, (n, b))) 
       push!(V, t1) 

       push!(I, index(ltc, (n, b))) 
       push!(J, index(ltc, (n, c))) 
       push!(V, t2)

       push!(I, index(ltc, (n, c))) 
       push!(J, index(ltc, (n1, b))) 
       push!(V, t2) 

    end
    H = sparse(I, J, V, ltc.N*ltc.U, ltc.N*ltc.U)
    return H + H'
end

