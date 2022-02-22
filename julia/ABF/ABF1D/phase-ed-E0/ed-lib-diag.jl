using LinearAlgebra, SparseArrays, KrylovKit
using LinearMaps
using Random
using DataFrames, CSV
# Custom modules
using ABF
using Lattice
using PN
include("./ed-params.jl") # read parameters from configuration file
include("./ed-lib-box1D.jl")

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


function generate_col_names(p::Params)
    return ["l"*string(li)*"_q"*string(qi) for li in p.l, qi in p.q]
end

function abf3d_scan(p::Params)
    col_str = generate_col_names(p)
    nt = Threads.nthreads()
    rng = [MersenneTwister(p.seed + i) for i in 1:nt]
    ltc = Lattice1D(p.L, 2)
    ltc_p = Lattice1D(p.L,  1)
    boxidx = [box_inds(ltc_p, li) for li in p.l]
    for j in p.th_ind_range[1]:p.th_ind_range[2]
        fn = "L$(p.L)_Th$(j)" #File name
        H, U = ham_fe(ltc, -2, 0, p.θ[j]) # Fully entangled hamiltonian
        H = convert.(ComplexF64, H)
        for jj in p.W_ind_range[1]:p.W_ind_range[2]
            df = [DataFrame(E = Float64[], r = Int64[]) for i in 1:nt]
            for t in 1:nt, k in 1:length(p.l), kk in 1:length(p.q)
                insertcols!(df[t], col_str[k, kk] => Float64[])
            end
            @Threads.threads for r in 1:p.R# Realizations
                er = true
                er_num = 0
                x = Threads.threadid()
                D_onsite = Diagonal(rand(rng[x], 2p.L) .- 0.5)
                D_phase = phase_dis(H, V = 1, rng = rng[x])
                @views H_prj = project(U'*(p.W[jj]*D_onsite + p.V*D_phase)*U)
                e, psi = eigen(Hermitian(Matrix(H_prj)))
                @views df_temp = DataFrame(E = round.(e[p.L÷2 + 1], sigdigits = 10), r = [r])
                #Push PN
                for k in 1:length(p.l), kk in 1:length(p.q)
                    @views insertcols!(df_temp, col_str[k, kk] => round.(compute_box_iprs(ltc_p, psi[:, p.L÷2 + 1], boxidx[k], q = p.q[kk]), sigdigits = 10))
                end
                append!(df[x], df_temp)
            end
            CSV.write(fn*"_W1_E1_purephase.csv", vcat(df...))
        end
    end
end
