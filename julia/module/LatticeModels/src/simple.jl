using Lattices
using SparseArrays

export chain, square, cubic

function chain(;N::Integer=10, t::F=1., periodic::Bool = true) where F
    ltc = Lattice1D(N, 1)
    I = Int64[]; J=Int64[]; V=F[]
    for n in 1:ltc.N    
       push!(I, index(ltc, (n, 1))) 
       push!(J, index(ltc, (n-1, 1))) 
       push!(V, -t) 

       push!(I, index(ltc, (n+1, 1))) 
       push!(J, index(ltc, (n, 1))) 
       push!(V, -conj(t)) 
    end
    return sparse(I, J, V, ltc.N, ltc.N)
end

function square(;M::Integer=10, N::Integer=10, t::F=1.) where F
    ltc = Lattice2D(M, N, 1)
    I = Int64[]; J=Int64[]; V=F[]
    for m in 1:ltc.M, n in 1:ltc.N    
       push!(I, index(ltc, (m,n,1))) 
       push!(J, index(ltc, (m,n-1,1))) 
       push!(V, -t) 

       push!(I, index(ltc, (m,n,1))) 
       push!(J, index(ltc, (m,n+1,1))) 
       push!(V, -conj(t)) 

       push!(I, index(ltc, (m,n,1))) 
       push!(J, index(ltc, (m-1,n,1))) 
       push!(V, -t) 

       push!(I, index(ltc, (m,n,1))) 
       push!(J, index(ltc, (m+1,n,1))) 
       push!(V, -conj(t)) 
    end
    return sparse(I, J, V, ltc.N*ltc.M, ltc.M*ltc.N)
end

function cubic(;L::Integer=10, M::Integer=10, N::Integer=10, t::F=1.) where F
    ltc = Lattice3D(L, M, N, 1)
    I = Int64[]; J=Int64[]; V=F[]
    for l in 1:ltc.L, m in 1:ltc.M, n in 1:ltc.N    
       push!(I, index(ltc, (l,m,n,1))) 
       push!(J, index(ltc, (l,m,n-1,1))) 
       push!(V, -t) 

       push!(I, index(ltc, (l,m,n,1))) 
       push!(J, index(ltc, (l,m,n+1,1))) 
       push!(V, -conj(t)) 

       push!(I, index(ltc, (l,m,n,1))) 
       push!(J, index(ltc, (l,m-1,n,1))) 
       push!(V, -t) 

       push!(I, index(ltc, (l,m,n,1))) 
       push!(J, index(ltc, (l,m+1,n,1))) 
       push!(V, -conj(t)) 

       push!(I, index(ltc, (l,m,n,1))) 
       push!(J, index(ltc, (l-1,m,n,1))) 
       push!(V, -t) 

       push!(I, index(ltc, (l,m,n,1))) 
       push!(J, index(ltc, (l+1,m,n,1))) 
       push!(V, -conj(t)) 
    end
    return sparse(I, J, V, ltc.L*ltc.N*ltc.M, ltc.L*ltc.M*ltc.N)
end
