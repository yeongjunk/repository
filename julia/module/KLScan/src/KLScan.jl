module KLScan

export shift_invert_linear_map, scan_kl
import Statistics: mean, std

using LinearAlgebra, SparseArrays, KrylovKit
using LinearMaps
using Random
using Lattices
using KLDivergence


"""
Linear map for A-IE
"""
function shift_invert_linear_map(A, E; c = 1., isherm = true)
    N = size(A, 1)
    F = lu(c*(A - E*I(N)))
    LinearMap{eltype(A)}((y, x) -> ldiv!(y, F, x), N, ismutating = true, ishermitian = isherm)
end

"""
Create Hamiltonian H = f(p; rng=GLOBAL_RNG) and compute nev number of eigenstates at the target energy E_c over the energy window E_del. 
R is number of realizations.
The box-counted PN is computed and averaged over realizations and the given energy window. 
From this, tau is computed
Optionally the parameter c is specified if an error occurs
"""

function scan_kl(f::Function, params, E_c, E_del; c=1., seed::Int = 1234, isherm::Bool = true, R::Int = 10, nev::Int= 10)
    L = params.L
    nt = Threads.nthreads()
    rng = [MersenneTwister(seed+i) for i in 1:nt]

    kls = [Float64[] for i in 1:nt]
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
                lmap = shift_invert_linear_map(H, E_c, c = c, isherm=isherm) 

                E_inv, psi, _ = eigsolve(lmap, n, nev, :LM, ishermitian=isherm, krylovdim = max(30, 2nev + 1));
                E = 1 ./ (c*real.(E_inv)) .+ E_c
             
                #---- Crop energies outside the energy bins ----#
                idx = findall(x -> (E_c-E_del) <= x <= (E_c+E_del), E)
                for i in 1:length(idx)
                    push!(E_full[x], E[idx[i]])
                end
                psi = reduce(hcat, psi[idx])
                #---- Compute KL1 ---#
                kls_temp = KL1(psi)

                # Push KL1 
                for i in 1:length(idx)
                    push!(E_full[x], E[idx[i]])
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

    E_full = vcat(E_full...)
    kls   = vcat(kls...) 

    kls_mean = mean(kls)
    E_mean = mean(E_full)

    return E_mean, kls_mean
end

end # module
