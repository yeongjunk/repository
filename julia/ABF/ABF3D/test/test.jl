using Lattice
using ABF
using LinearAlgebra
using SparseArrays
using Plots
using Arpack
using Random
using LinearMaps
using LogExpFunctions
using Statistics
using KrylovKit
gr(fmt=:png)

function phase_dis(H; V = 1., rng = nothing)
    D = convert.(ComplexF64, H)
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
              vals[i] = im*V*vals[i]*(rand(rng) .- 0.5)
          elseif row <= j
              vals[i] = 0.
          end
       end
    end
    return D + D'
end

function construct_linear_map(A)
    F = factorize(A)
    LinearMap{eltype(A)}((y, x) -> ldiv!(y, F, x), size(A, 1), ismutating = true)
end

function box_idx2d(ltc, b)
    box_inds = Array{Int64, 1}[]
    L = ltc.M
    for y in 1:b:(L-b+1), z in 1:b:(L-b+1)
        idx = Int64[]
        for m in 0:b-1, n in 0:b-1
            push!(idx, index(ltc, (y + m, z + n, 1)))
        end
        push!(box_inds, idx)
    end
    return box_inds
end

L = 100
ltc = Lattice2D(L, L, 2)
H, U = ham_fe(ltc, -2, 0, 0.25);
D = phase_dis(H)
H_sf = project(U'*(H + D)*U);
H_sf = sparse(Hermitian(H_sf))
droptol!(H_sf, 1E-12);
@time e_inv, psi, info = eigsolve(construct_linear_map(H_sf .- 0.001I(size(H_sf, 1))), size(H_sf, 1), 10, :LM, issymmetric = true, krylovdim = max(30, 21));
nnz(H_sf)

L = 22
ltc_3d = Lattice3D(L, L, L, 2)
H, U = ham_fe(ltc_3d, -2, 0, 0.25);
D = Diagonal(rand(size(H, 1)) .- 0.5)
H_sf = project(U'*(H + D)*U);
droptol!(H_sf, 1E-12);
H_sf = Hermitian(H_sf)
@time e_inv, psi, info = eigsolve(construct_linear_map(H_sf), size(H_sf, 1), 10, :LM, issymmetric = true, krylovdim = max(30, 21));


nnz(H_sf)
