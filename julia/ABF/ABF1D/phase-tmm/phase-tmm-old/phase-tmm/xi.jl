using Random
using LinearAlgebra
using SparseArrays
using Lattice, ABF
using Parameters
include("./params.jl")
@doc """
Unlike off_phase_dis2!, this returns 'phase disorder matrix' only, but
it does not add to the original hamiltonian.
"""
function off_phase_dis2(H, r_off)
    H_copy = copy(H)
    rows = rowvals(H_copy)
    vals = nonzeros(H_copy)

    for j = 1:size(H_copy, 1)
       for i in nzrange(H_copy, j)
          if j > rows[i]
              vals[i] = vals[i]*im*r_off[i]
              H_copy[j,rows[i]] = conj(vals[i])
          elseif j == rows[i]
              vals[i] = 0.
          end
      end
    end
    dropzeros!(H_copy)
    return H_copy
end


@doc """
transfer matrix of nth supercell...
"""
function tm(E, t_1, t_2, t_m1, t_m2, tt_1, tt_2, tt_m1, tt_m2)
    TM = Array{ComplexF64}(undef, 4, 4)
    H_n_0 = [E -t_m1; -tt_1 E]
    T_n_1 = [t_2 t_1; 0 tt_2]
    T_n_m1 = [t_m2 0; tt_m1 tt_m2]
    TM[1:2, 1:2] = inv(T_n_1)*H_n_0
    TM[1:2, 3:4] = -inv(T_n_1)*T_n_m1
    TM[3:4, 1:2] = I(2)
    TM[3:4, 3:4] .= 0
    return TM
end

function ham_sf_ph(H_fe, D_dis, T_dis, U, W, V)
    H_dis = H_fe + V*T_dis
    H_dis = U'*H_dis*U
    H_out = projection(H_dis)
    return H_out
end

@doc """
Add onsite & phase disorder -> detangle -> project(scale free, normalize bandwidth)
"""
function ham_sf_ph_on(H_fe, D_ph, D_on, U, W, V)
    H_dis = H_fe + V*D_ph + W*D_on
    H_dis = U'*H_dis*U
    H_out = projection(H_dis)
    return H_out
end

function correlation(x)
    cor = Float64[]
    for i in (length(x)÷4+1):(length(x)÷4*3)
        push!(cor, 0.)
    end
    for i in (length(x)÷4+1):(length(x)÷4*3)
        for r in 1:(length(x)÷4)
            cor[r] += abs(x[i+r] - x[i])
        end
    end
    return cor
end

function loc_length(;E = 0., V = 1., W = 0., θ = 0.25, δ = 0., N = 10000, q = 4, seed = 1234)
    ltc = Lattice1D(N, 2)
    H_fe, U = ham_fe(ltc, -1., 1., θ);
    H_fe = convert.(ComplexF64, H_fe)
    rng = MersenneTwister(seed)
    r_on = rand(rng, 2N) .- 0.5
    r_off = rand(rng, 12N) .- 0.5

    D = spdiagm(0=>r_on)
    T = off_phase_dis2(H_fe, r_off);
    H_sf = Hermitian(project(U'*(V*T + W*D)*U));

    # Wavefunction with arbitrary precision
    ψⱼ = convert.(ComplexF64, I(4))
    ΣlnR = 0.;


    for j in 2:2:N-2
        t₁ⱼm1 = H_sf[mod1(j-1, N),mod1(j, N)]
        t₂ⱼm1 = H_sf[mod1(j-1, N),mod1(j+1, N)]
        tm₁ⱼm1 = H_sf[mod1(j-1, N),mod1(j-2, N)]
        tm₂ⱼm1 = H_sf[mod1(j-1, N),mod1(j-3, N)]

        t₁ⱼ = H_sf[mod1(j,N),mod1(j+1, N)]
        t₂ⱼ = H_sf[mod1(j,N),mod1(j+2, N)]
        tm₁ⱼ = H_sf[mod1(j,N),mod1(j-1, N)]
        tm₂ⱼ = H_sf[mod1(j,N),mod1(j-2, N)]

        eⱼ = H_sf[mod1(j,N),mod1(j,N)]

        Tⱼ = tm(E, t₁ⱼ, t₂ⱼ, tm₁ⱼ, tm₂ⱼ, t₁ⱼm1, t₂ⱼm1,  tm₁ⱼm1, tm₂ⱼm1)
        ψⱼ = Tⱼ*ψⱼ

        # QR
        if j%q == 0
            ψⱼ_new, R = qr!(ψⱼ)
            ψⱼ .= Matrix(ψⱼ_new)
            ΣlnR += log(abs.(R[2, 2]))
            ψⱼ .= Matrix{ComplexF64}(ψⱼ)
        end
    end
    return Float64(2N/ΣlnR)
end

function loc_length(p::Params)
    @unpack E, θ, V, W, q, N, seed = p
    return loc_length(E = E, θ = θ, V = V, W = W, q = q, N = N, seed = seed)
end
