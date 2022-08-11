using Lattice
using ABF
using Random
using SparseArrays
using LinearAlgebra  
using Distributions

function hop_dis(H; V = 1., rng = nothing)
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
              vals[i] = V*vals[i]*(rand(rng, F) .- 0.5)
          elseif row <= j
              vals[i] = F(0.)
          end 
       end 
    end 
    return D + D'
end

function square_pbc(ltc::Lattice2D)
    I = Int64[]; J = Int64[];
    V = Float64[]; 
    for m in 1:ltc.M, n in 1:ltc.N
        push!(I, index(ltc, (m,n,1)))
        push!(J, index(ltc, (m,n+1,1)))
        push!(V, -1.)

        push!(I, index(ltc, (m,n,1)))
        push!(J, index(ltc, (m+1,n,1)))
        push!(V, -1.)

        push!(I, index(ltc, (m,n,1)))
        push!(J, index(ltc, (m-1,n,1)))
        push!(V, -1.)

        push!(I, index(ltc, (m, n,1)))
        push!(J, index(ltc, (m, n-1,1)))
        push!(V, -1.)
    end
    return sparse(I, J, V, ltc.M*ltc.N, ltc.M*ltc.N)
end

function square_obc(ltc::Lattice2D)
    I = Int64[]; J = Int64[];
    V = Float64[]; 
    for m in 1:ltc.M, n in 1:ltc.N
        if n < ltc.N
            push!(I, index(ltc, (m,n,1)))
            push!(J, index(ltc, (m,n+1,1)))
            push!(V, -1.)
        end
        if m < ltc.N
            push!(I, index(ltc, (m,n,1)))
            push!(J, index(ltc, (m+1,n,1)))
            push!(V, -1.)
        end

        if m > 1
            push!(I, index(ltc, (m,n,1)))
            push!(J, index(ltc, (m-1,n,1)))
            push!(V, -1.)
        end

        if n > 1
            push!(I, index(ltc, (m,n,1)))
            push!(J, index(ltc, (m,n-1,1)))
            push!(V, -1.)
        end
    end
    return sparse(I, J, V, ltc.M*ltc.N, ltc.M*ltc.N)
end
