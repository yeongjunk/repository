module KLScan

export shift_invert_linear_map, scan_kl
import Statistics: mean, std

using LinearAlgebra, SparseArrays, KrylovKit
using LinearMaps
using Random
using Lattices
using KLDivergence

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

function shift_invert_linear_map(A, E; c = 1.0, isherm = true)
    N = size(A, 1)
    F = factorize(c*(A - E*I(N)))
    LinearMap{eltype(A)}((y, x) -> ldiv2!(y, F, x), N, ismutating = true, ishermitian = isherm)
end

function scan_kl(f::Function, params, E_c, E_del; c = 1.0, rng::AbstractRNG = Random.GLOBAL_RNG, isherm::Bool = true, R::Int = 10, nev::Int= 10, kwarg_eig...)
    # Initalize thread safe rngs
    nt         = Threads.nthreads()
    masterseed = rand(rng, 100:99999999999)
    rngs       = generateParallelRngs(rng, nt)

    # Initialize output array
    kls        = [Float64[] for i in 1:nt]
    Es         = [Float64[] for i in 1:nt]
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
                lmap          = shift_invert_linear_map(H, E_c, c = c, isherm=isherm) 
                E_inv, psi, _ = eigsolve(lmap, n, nev, :LM; kwarg_eig...)
                E_temp        = 1 ./ (c*real.(E_inv)) .+ E_c 

                #---- Crop energies outside the energy bins ----#
                idx      = findall(x -> (E_c-E_del) <= x <= (E_c+E_del), E_temp)
                psi      = reduce(hcat, psi[idx])
                kls_temp = KL1(psi)

                for i in 1:length(idx)
                    push!(Es[x], E_temp[idx[i]])
                    if i != length(idx) 
                        push!(kls[x], kls_temp[i])
                    end
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

    Es       = vcat(Es...)
    Es_mean  = mean(Es)
    kls      = vcat(kls...) 
    kls_mean = mean(kls)
    return Es_mean, kls_mean
end

end # module
