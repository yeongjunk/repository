module TAU

export shift_invert_linear_map, compute_tau
import Statistics: mean, std

using LinearAlgebra, SparseArrays, KrylovKit
using LinearMaps
using Random
using DataFrames
using ArgParse, JSON
using Lattices
using PN
using GIPR
using Statistics

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
function compute_tau(f::Function, params, E_c, E_del, R::Int, nev::Int, ltc::Lattice, L::Int, l::Int; c=1., seed::Int = 1234, isherm = true)
    nt = Threads.nthreads()
    rng = [MersenneTwister(seed) for i in 1:nt]

    # Allocate Empty dataframe
    df = [DataFrame(E = Float64[], iprs = Float64[]) for i in 1:nt]

    box_inds = box_indices(ltc, l)
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
                psi = reduce(hcat, psi[idx])


                #---- Compute box-counted IPR ---#
                p = abs2.(psi)            
                p_coarse = box_coarse(p, box_inds)
                iprs = compute_iprs(p_coarse, density=true)
                df_temp = DataFrame(E = E[idx], iprs = iprs)
                # Push PN
                append!(df[x], df_temp)
                er = false
            catch e
                println("There was an error: ", e.msg)
                num_try += 1
                continue
            end # try
        end # while
    end # for loop over Realization

    df_full = DataFrame(E = Float64[], iprs = Float64[])
    for x in 1:nt
        append!(df_full, df[x])
    end

    E_mean = mean(df_full.E)
    ipr_mean = mean(df_full.iprs)
    
    tau = log(ipr_mean)/log(l/L)

    return E_mean, tau 
end

end # module
