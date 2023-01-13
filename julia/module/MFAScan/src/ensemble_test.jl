module EnsembleTest

export shift_invert_linear_map, scan_ταf, mt_scan_ταf, MFAParameters
import Statistics: mean, std

using Statistics
using LinearAlgebra, SparseArrays, KrylovKit
using LinearMaps
using Random
using Lattices
using PN
using GIPR
using MFA.Ensemble

dropmean(A; dims=:) = dropdims(mean(A; dims=dims); dims=dims)

"""
    generateParallelRngs(rng::AbstractRNG, n::Integer;reSeed=false)
For multi-threaded models, return n independent random number generators (one per thread) to be used in threaded computations.
Note that each ring is a _copy_ of the original random ring. This means that code that _use_ these RNGs will not change the original RNG state.
Use it with `rngs = generateParallelRngs(rng,Threads.nthreads())` to have a separate rng per thread.
By default the function doesn't re-seed the RNG, as you may want to have a loop index based re-seeding strategy rather than a threadid-based one (to guarantee the same result independently of the number of threads).
If you prefer, you can instead re-seed the RNG here (using the parameter `reSeed=true`), such that each thread has a different seed. Be aware however that the stream  of number generated will depend from the number of threads at run time.
"""
function generateParallelRngs(rng::AbstractRNG, n::Integer;reSeed=false)
    if reSeed
        seeds = [rand(rng,100:18446744073709551615) for i in 1:n] # some RNGs have issues with too small seed
        rngs  = [deepcopy(rng) for i in 1:n]
        return Random.seed!.(rngs,seeds)
    else
        return [deepcopy(rng) for i in 1:n]
    end
end

function ldiv2!(y, F, x)
    y .= F\x
end

"""
Linear map for A-IE
"""
function shift_invert_linear_map(A, E; c = 1., isherm = true)
    N = size(A, 1)
    F = factorize(c*(A - E*I(N)))
    LinearMap{eltype(A)}((y, x) -> ldiv2!(y, F, x), N, ismutating = true, ishermitian = isherm)
end

"""
Create Hamiltonian H = f(p; rng=GLOBAL_RNG) and compute nev number of eigenstates at the target energy E_c over the energy window E_del. 
R is number of realizations.
The box-counted PN is computed and averaged over realizations and the given energy window. 
From this, tau is computed
Optionally the parameter c is specified if an error occurs
"""
function scan_ταf(f::Function, params, ε::Float64, Δε::Float64, p_MFA::MFAParameters; nev = 10, rng = Random.GLOBAL_RNG, isherm::Bool = true, noaverage=false, R::Int = 10, kwargs_eig...)
    nt = Threads.nthreads() 
    masterseed = rand(rng, 100:99999999999)
    rngs       = generateParallelRngs(rng, nt)

    gipr  = [Array{Float64}[] for i in 1:nt]
    μqlnμ = [Array{Float64}[] for i in 1:nt]
    E     = [Float64[] for i in 1:nt]
    #----------------------------- DIAGONALIZATION -----------------------------#
    for r in 1:R # Realizations
        er = true
        num_try = 0
        while er
            num_try == 5 && error("5 attempts faild.")
            try
                x = Threads.threadid()
                tsrng = rngs[x]
                Random.seed!(tsrng,masterseed+r*40)
                #---- Create the model ----#
                H = f(params, rng = tsrng)
                n = size(H, 1)

                #---- Lanczos method with shift-and-invert method ----#
                lmap = shift_invert_linear_map(H, ε, isherm=isherm) 
                E_temp, psi, _ = eigsolve(lmap, n, nev, :LM; kwargs_eig...) 
                @. E_temp = 1 / real(E_temp) + ε 
                E_temp = convert.(Float64, E_temp)
                #---- Crop energies outside the energy bins ----#
                idx = findall(x -> (ε - Δε) <= x <= (ε + Δε), E_temp)

                #---- Compute GIPR ---#
                for i in 1:length(idx)
                    gipr_temp, μqlnμ_temp = compute_gipr_2(p_MFA, psi[idx[i]])
                    push!(E[x], E_temp[idx[i]])
                    push!(gipr[x], gipr_temp)
                    push!(μqlnμ[x], μqlnμ_temp)
                end

                er = false
            catch e
                println("There was an error: ")
                @error "ERROR: " exception=(e, catch_backtrace())
                num_try += 1
                continue
            end # try
        end # while
    end # for loop over Realization
    E = vcat(E...)
    gipr   = vcat(gipr...) 
    gipr   = reshape.(gipr, 1, size(gipr[1])...)
    gipr   = vcat(gipr...) 

    μqlnμ  = vcat(μqlnμ...) 
    μqlnμ  = reshape.(μqlnμ, 1, size(μqlnμ[1])...)
    μqlnμ  = vcat(μqlnμ...) 

    E_mean = mean(E)
    gipr_mean = dropmean(gipr, dims = 1)
    μqlnμ_mean = dropmean(μqlnμ, dims = 1)
      
    τ, α, f_α = compute_ταf(p_MFA, gipr, μqlnμ)

    if !noaverage
        return E_mean, gipr_mean, μqlnμ_mean, τ, α, f_α 
    else
        return E, gipr, μqlnμ, τ, α, f_α
    end
    
end



"""
Create Hamiltonian H = f(p; rng=GLOBAL_RNG) and compute nev number of eigenstates at the target energy E_c over the energy window E_del. 
R is number of realizations.
The box-counted PN is computed and averaged over realizations and the given energy window. 
From this, tau is computed
Optionally the parameter c is specified if an error occurs
"""
function mt_scan_ταf(f::Function, params, ε::Float64, Δε::Float64, p_MFA::MFAParameters; nev = 10, rng = Random.GLOBAL_RNG, isherm::Bool = true, noaverage=false, R::Int = 10, kwargs_eig...)
    nt = Threads.nthreads() 
    masterseed = rand(rng, 100:99999999999)
    rngs       = generateParallelRngs(rng, nt)

    gipr  = [Array{Float64}[] for i in 1:nt]
    μqlnμ = [Array{Float64}[] for i in 1:nt]
    E     = [Float64[] for i in 1:nt]
    #----------------------------- DIAGONALIZATION -----------------------------#
    @Threads.threads for r in 1:R # Realizations
        er = true
        num_try = 0
        while er
            num_try == 5 && error("5 attempts faild.")
            try
                x = Threads.threadid()
                tsrng = rngs[x]
                Random.seed!(tsrng,masterseed+r*40)
                #---- Create the model ----#
                H = f(params, rng = tsrng)
                n = size(H, 1)

                #---- Lanczos method with shift-and-invert method ----#
                lmap = shift_invert_linear_map(H, ε, isherm=isherm) 
                E_temp, psi, _ = eigsolve(lmap, n, nev, :LM; kwargs_eig...) 
                @. E_temp = 1 / real(E_temp) + ε 
                E_temp = convert.(Float64, E_temp)
                #---- Crop energies outside the energy bins ----#
                idx = findall(x -> (ε - Δε) <= x <= (ε + Δε), E_temp)

                #---- Compute GIPR ---#
                for i in 1:length(idx)
                    gipr_temp, μqlnμ_temp = compute_gipr_2(p_MFA, psi[idx[i]])
                    push!(E[x], E_temp[idx[i]])
                    push!(gipr[x], gipr_temp)
                    push!(μqlnμ[x], μqlnμ_temp)
                end

                er = false
            catch e
                println("There was an error: ")
                @error "ERROR: " exception=(e, catch_backtrace())
                num_try += 1
                continue
            end # try
        end # while
    end # for loop over Realization
    E = vcat(E...)
    gipr   = vcat(gipr...) 
    gipr   = reshape.(gipr, 1, size(gipr[1])...)
    gipr   = vcat(gipr...) 

    μqlnμ  = vcat(μqlnμ...) 
    μqlnμ  = reshape.(μqlnμ, 1, size(μqlnμ[1])...)
    μqlnμ  = vcat(μqlnμ...) 

    E_mean = mean(E)
    gipr_mean = dropmean(gipr, dims = 1)
    μqlnμ_mean = dropmean(μqlnμ, dims = 1)
      
    τ, α, f_α = compute_ταf(p_MFA, gipr, μqlnμ)

    if !noaverage
        return E_mean, gipr_mean, μqlnμ_mean, τ, α, f_α 
    else
        return E, gipr, μqlnμ, τ, α, f_α
    end
    
end
end # module
