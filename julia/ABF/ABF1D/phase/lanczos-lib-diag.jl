using LinearAlgebra, SparseArrays, KrylovKit
using LinearMaps
using Random
using DataFrames, CSV
# Custom modules
using ABF
using Lattice
using PN
include("./lanczos-params.jl") # read parameters from configuration file
include("./lanczos-lib-box1D.jl")

# Functions related to lanczos matrix diagonalization
function construct_linear_map(A)
    F = factorize(A)
    LinearMap{eltype(A)}((y, x) -> ldiv2!(y, F, x), size(A, 1), ismutating = true, ishermitian = true)
end
function ldiv2!(y, F, x)
    y .= F\x
end

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

function estimate_bw(p, θ, V, W, L, rng)
    ltc = Lattice1D(L, 2)
    H, U = ham_fe(ltc, -2, 0, θ) # Fully entangled hamiltonian
    H = convert.(ComplexF64, H)
    D_phase = phase_dis(H, V = 1, rng = rng)
    D_onsite = Diagonal(rand(rng, 2L) .- 0.5)

    @views H_prj = project(U'*(V*D_phase +W*D_onsite)*U)
    vals, psi, info = eigsolve(Hermitian(H_prj), size(H_prj, 1), 1, :LM, ishermitian = true, krylovdim = 30)
    return 2*abs(maximum(vals))
end

# @doc"""
# p: Params
# θ: angle
# W: disorder strength
# L: system size for estimating bandwidth.
# E_crop: crop the end of the spectrum. Recommend 0.9
# E_bin_width: portion given by 1/E_bin_width of the small energy window around energy center.
# rng: random number generator
# """
function auto_energy_params(p, θ,V, W, L, E_crop, E_bin_width, rng)
    BW = estimate_bw(p, θ, V, W, L, rng)
    E_c = range(0.0001, BW/2*E_crop, length = length(p.E))
    E_del = (E_c[2] - E_c[1])/E_bin_width
    return BW, E_c, E_del
end

function energy_param_generator(p, θ, V, W, L, E_crop, E_bin_width, rng)
    if p.bw_auto
        BW, E_c, E_del = auto_energy_params(p, θ, V, W, L, E_crop, E_bin_width, rng)
    else
        if length(p.E) != 1
            BW = p.E[end] - p.E[1]
            E_c = p.E
            E_del = (E_c[2] - E_c[1])/p.E_bin_width
        elseif length(p.E) == 1
            BW = estimate_bw(p, θ, V, W, L, rng)
            E_c = p.E[1] + 1E-12
            E_del = BW/p.E_bin_width
        end
    end
    return BW, E_c, E_del
end

function print_info(nt, BW, E_c, E_del, p)
    println("Number of BLAS threads: ",LinearAlgebra.BLAS.get_num_threads())
    println("Number of threads: $(nt)")
    println("Realization / nev: ", p.R/p.nev)
    println("Auto bandwidth settings: ", p.bw_auto)
    println("Bandwidth estimated: $(BW)")
    println("Energy centers: ", E_c)
    println("Energy bin width: $(E_del)")
    println("Thetas: ", p.θ)
end

function generate_col_names(p::Params)
    return ["l"*string(li)*"_q"*string(qi) for li in p.l, qi in p.q]
end

@time function abf3d_scan(p::Params)
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
            BW, E_c, E_del = energy_param_generator(p, p.θ[j], p.V, p.W[jj], 5000, 0.9, 8, rng[1])
            print_info(nt, BW, E_c, E_del, p)
            for jjj in p.E_ind_range[1]:p.E_ind_range[2]
                df = [DataFrame(E = Float64[], r = Int64[]) for i in 1:nt]
                for t in 1:nt, k in 1:length(p.l), kk in 1:length(p.q)
                    insertcols!(df[t], col_str[k, kk] => Float64[])
                end
                @Threads.threads for r in 1:p.R÷p.nev# Realizations
                    er = true
                    er_num = 0
                    while er
                        er_num == 5 && error("5 attempts faild.")
                        try
                            x = Threads.threadid()
                            D_onsite = Diagonal(rand(rng[x], 2p.L) .- 0.5)
                            D_phase = phase_dis(H, V = 1, rng = rng[x])

                            @views H_prj = project(U'*(p.W[jj]*D_onsite + p.V*D_phase)*U)
                            droptol!(H_prj, 1E-13)
                            e_inv, psi, info = eigsolve(construct_linear_map(Hermitian(p.L*(H_prj - E_c[jjj]*I(size(H_prj, 1))))), size(H_prj, 1),
                                p.nev, :LM, ishermitian = true, krylovdim = max(30, 2p.nev + 1));
                            e = 1 ./ (p.L*real.(e_inv)) .+ E_c[jjj]
                            psi = reduce(hcat, psi)

                            #Crop energies that are outside the energy bins
                            idx = findall(x -> (E_c[jjj] - E_del) <  x < (E_c[jjj] + E_del), e)
                            @views df_temp = DataFrame(E = round.(e[idx], sigdigits = 10), r = fill(r, length(idx)))
                            #Push PN
                            for k in 1:length(p.l), kk in 1:length(p.q)
                                @views insertcols!(df_temp, col_str[k, kk] => round.(compute_box_iprs(ltc_p, psi[:, idx], boxidx[k], q = p.q[kk]), sigdigits = 10))
                            end
                            append!(df[x], df_temp)
                            er = false
                        catch e
                            println("There was an error: ", e.msg)
                            er_num += 1
                            continue
                        end
                    end
                end
                CSV.write(fn*"_W$(lpad(jj, 2, "0"))_E$(lpad(jjj, 2, "0")).csv", vcat(df...))
            end
        end
    end
end
