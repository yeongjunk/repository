module Typical 

export shift_invert_linear_map, scan_ταf
import Statistics: mean, std

using Statistics
using LinearAlgebra, SparseArrays, KrylovKit
using LinearMaps
using Random
using DataFrames
using ArgParse, JSON
using Lattices
using PN
using GIPR
using MFA.Typical

dropmean(A; dims=:) = dropdims(mean(A; dims=dims); dims=dims)
"""
Linear map for A-IE
"""
function shift_invert_linear_map(A, E; c = 1., isherm = true)
    N = size(A, 1)
    F = lu(c*(A - E*I(N)))
    LinearMap{eltype(A)}((y, x) -> ldiv!(y, F, x), size(A, 1), ismutating = true, ishermitian = isherm)
end
"""
compute_tau(f, p, E_c, E_del, R, nev, L, l; c = 1., seed = 1234) 
Create Hamiltonian H = f(p; rng=GLOBAL_RNG) and compute nev number of eigenstates at the target energy E_c over the energy window E_del. 
R is number of realizations.
The box-counted PN is computed and averaged over realizations and the given energy window. 
From this, tau is computed
Optionally the parameter c is specified if an error occurs
"""

function scan_ταf(f::Function, params, E_c, E_del, ltc::Lattice; c=1., seed::Int = 1234, isherm::Bool = true, l::Vector{Int} = [1],  q::Vector{Float64} = [2.], R::Int = 10, nev::Int= 10)
    L = ltc.N
    nt = Threads.nthreads()
    rng = [MersenneTwister(seed) for i in 1:nt]

    p_MFA = MFAParameters(ltc, l, q)
    prepare_MFA!(p_MFA)
    τ = [Array{Float64}[] for i in 1:nt]
    α = [Array{Float64}[] for i in 1:nt]
    f_α = [Array{Float64}[] for i in 1:nt]
    E_full = [Float64[] for i in 1:nt]
    #----------------------------- DIAGONALIZATION -----------------------------#
    @Threads.threads for r in 1:R # Realizations
        er = true
        num_try = 0
        while er
            num_try == 5 && error("5 attempts faild.")
            try
                x = Threads.threadid()
                #---- Create the model ----#
                H = f(params, rng = rng[x])
                n = size(H, 1)

                #---- Lanczos method with shift-and-invert method ----#
                lmap = shift_invert_linear_map(H, E_c, c = c) 

                E_inv, psi, _ = eigsolve(lmap, n, nev, :LM, ishermitian=isherm, krylovdim = max(30, 2nev + 1));
                E = 1 ./ (c*real.(E_inv)) .+ E_c
             
                #---- Crop energies outside the energy bins ----#
                idx = findall(x -> (E_c-E_del) <= x <= (E_c+E_del), E)
                for i in 1:length(idx)
                    push!(E_full[x], E[idx[i]])
                end
                psi = reduce(hcat, psi[idx])

                #---- Compute GIPR ---#
                τ_temp = Array{Float64}(undef, size(psi, 2), length(p_MFA.q), length(p_MFA.l))  
                α_temp = similar(τ_temp) 
                f_α_temp = similar(τ_temp) 

                for i in 1:size(psi, 2)
                    τ_temp[i, :, :], α_temp[i, :, :], f_α_temp[i, :, :] = compute_ταf(p_MFA, psi[:, i])
                end

                # Push GIPR 
                push!(τ[x], τ_temp)
                push!(α[x], α_temp)
                push!(f_α[x], f_α_temp)
                er = false
            catch e
                println("There was an error: ", e.msg)
                num_try += 1
                continue
            end # try
        end # while
    end # for loop over Realization
    E_full = reduce(vcat, E_full)
    E_mean = mean(E_full)

    τ = reduce(vcat, τ) 
    τ = reduce(vcat, τ) 
    α = reduce(vcat, α) 
    α = reduce(vcat, α) 
    f_α = reduce(vcat, f_α) 
    f_α = reduce(vcat, f_α) 

    return E_mean, dropmean(τ, dims = 1), dropmean(α, dims = 1), dropmean(f_α, dims = 1)
end

end # module
