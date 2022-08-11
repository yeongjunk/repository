using Lattice
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


function iter_tuple(i::Integer, j::Integer, minval::Integer, X::Integer, Y::Integer)
    if j != Y
        return (i, j+1)
    elseif i != X
        return (i + 1, minval)
    elseif i == X && j == Y
        return (Inf, Inf)
    end
end

"""
check the range of hoppings
"""
function hop_checker(ltc, H)
    if !(ltc.M >= 20 && ltc.N >= 20)
        error("The lattice size is not big enough")
    end
    U = ltc.U
    M = ltc.M ÷ 2
    N = ltc.N ÷ 2
    
    idx1 = [index(ltc, (M, N, u)) for u in 1:U]
    flag = true
    i = -10
    j = -10
    nzs = Tuple{Int64, Int64}[];
    while i != Inf && j != Inf
        idx2 = [index(ltc, (M+i, N+j, u)) for u in 1:U]
        iszeroblock = all(abs.(H[idx1, idx2]) .< 1E-13)
        if !iszeroblock
            push!(nzs, (i, j))
        end
        i, j = iter_tuple(i, j, -10, ltc.M - M, ltc.N - N)
    end
    return nzs
end  


"""
check the range of hoppings
"""
function pbc_remover!(ltc, H)
    hops = hop_checker(ltc, H)
    for i in 1:length(hops)
        x, y = hops[i]
        for m in 1:ltc.M, n in 1:ltc.N
            if (m + x > ltc.M) || (n + y > ltc.N) || (m + x < 1) || (n + y < 1)    
                idx1 = [index(ltc, (m, n, u)) for u in 1:ltc.U]
                idx2 = [index(ltc, (m+x,n+y, u)) for u in 1:ltc.U]
                H[idx1, idx2] = zeros(eltype(H), ltc.U, ltc.U)
                H[idx2, idx1] = zeros(eltype(H), ltc.U, ltc.U)
            end
        end
    end
    dropzeros!(H)
end  


function ham_sf_obc(;L::Integer = 10, V = 1, θ::F = 0.25,rng = Random.GLOBAL_RNG) where F <: AbstractFloat
    ltc = Lattice2D(L, L, 2)
    ltc_p = Lattice2D(L, L, 1)
    H, U = ham_fe(ltc, F(-2.), F(0.), θ)
    D = phase_dis(H, V = V, rng = rng)
    H_sf = project(U'*D*U)
    pbc_remover!(ltc_p,H_sf)
    return H_sf
end

function ham_sf_pbc(;L::Integer = 10, V = 1, θ::F = 0.25,rng = Random.GLOBAL_RNG) where F <: AbstractFloat
    ltc = Lattice2D(L, L, 2)
    H, U = ham_fe(ltc, F(-2.), F(0.), θ)
    D = phase_dis(H, V = V, rng = rng)
    H_sf = project(U'*D*U)
    return H_sf
end

function compute_psi(; L::Integer = 1001,θ::F = 0.25, rng::AbstractRNG = Random.GLOBAL_RNG) where F <: AbstractFloat
    H = ham_sf_obc(L = L, θ = θ, rng = rng)
    H = Complex{F}.(BandedMatrix(H, (2L+2, 2L+2)));
    Y = Array(@view H[:, 1])
    A = H[:, 2:end]
    X = A\Y
    pushfirst!(X, 1.)    
    return X
end
