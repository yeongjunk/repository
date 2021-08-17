using ABF
using Lattice
using PN
using LinearAlgebra
using SparseArrays
using Random
using Plots
using Arpack
using Statistics
using DataFrames

using SparseArrays
using LinearAlgebra
using Random
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

function hop_dis(H; V = 1., rng = nothing)
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
              vals[i] = V*vals[i]*(rand(rng) .- 0.5)
          elseif row <= j
              vals[i] = 0.
          end
       end
    end
    return D + D'
end

function GIPR(psi; q = 2)
    return compute_iprs(psi, q = q)
end

function compute_fα(p, b, L; q = 2.)
    μ = p.^q
    Z_q = sum(μ)
    for i in 1:size(p, 2)
        μ[:, i] ./= Z_q[i]
    end
    mu_dummy = replace(μ, 0. => 1.)
    return (1/log(b/L))*sum(μ.*log.(mu_dummy))
end

function compute_α(p, b, L; q = 2.)
    μ = p.^q
    Z_q = sum(μ)
    for i in 1:size(p, 2)
        μ[:, i] ./= Z_q[i]
    end
    return (1/log(b/L))*sum(μ.*log.(abs.(p)))
end



function box_idx3d(ltc, b)
    box_inds = Array{Int64, 1}[]
    L = ltc.L
    for x in 1:b:(L-b+1), y in 1:b:(L-b+1), z in 1:b:(L-b+1)
        idx = Int64[]
        for l in 0:b-1, m in 0:b-1, n in 0:b-1
            push!(idx, index(ltc, (x + l, y + m, z + n, 1)))
        end
        push!(box_inds, idx)
    end
    return box_inds
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

L_full = 20:20:200

df = DataFrame(q = Float64[], a = Float64[], fa = Float64[], tau = Float64[], L = Float64[])
for (k, L) in enumerate(L_full)
    bs = [1;]
    qs = 0:0.1:3.

    ltc = Lattice2D(L, L, 2)
    l_p = Lattice2D(L, L, 1)
    H, U = ham_fe(ltc, -1, 1, 0.25);
    # D = spdiagm(0 => (rand(size(H, 1)) .- 0.2));
    D = phase_dis(H)
    E, psi = eigen(Hermitian((Matrix(project(U'*(H + D)*U)))), 0.99, 1.01);
    R = size(psi, 2)
    α = Array{Float64}(undef, length(qs), length(bs), R)
    fα = similar(α, length(qs), length(bs), R)
    for r in 1:R
        for i in 1:length(bs)
            inds = box_idx2d(l_p, bs[i])
            p = [sum(psi[inds[i], r]) for i in 1:length(inds)]
            p = abs.(psi[:, r]).^2
            for j in 1:length(qs)
                α[j, i, r] = compute_α(p, bs[i], L, q = qs[j])[1, 1]
                fα[j, i, r] = compute_fα(p, bs[i], L, q = qs[j])[1, 1]
            end
        end
    end

    α = α[:,1,:];
    fα = fα[:,1,:];

    α = mean(α, dims = 2);
    fα = mean(fα, dims = 2);
    τ_q = α.*qs .- fα;
    L_rep = repeat([L], length(qs))

    df_temp = DataFrame(q = vec(qs), a = vec(α), fa = vec(fα), tau = vec(τ_q), L = L_rep)
    append!(df, df_temp)
end
