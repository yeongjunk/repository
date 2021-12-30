using Random
using LinearAlgebra
using SparseArrays
using StaticArrays
@doc """
Hamiltonian of isolated strip used for RGF.
"""
function isolated_strip!(H_n; M::Integer = 3, W::Real = 1., t::Real = 1., rng = Random.GLOBAL_RNG)
    @assert M >= 3 # Required due to periodic boundary condition
    for i in 1:M
        H_n[i, mod1(i+1, M)] = t
        H_n[i, mod1(i-1, M)] = t
        H_n[mod1(i+1, M), i] = t
        H_n[mod1(i-1, M), i] = t
    end
end

@doc """
Hamiltonian of isolated strip used for RGF.
"""
function isolated_strip(; M::Integer = 3, W::Real = 1., t::Real = 1., rng = Random.GLOBAL_RNG)
    @assert M >= 3 # Required due to periodic boundary condition
    H_n = spzeros(Float64, M^2, M^2)
    for i in 1:M
        H_n[i, mod1(i+1, M)] = t
        H_n[i, mod1(i-1, M)] = t
        H_n[mod1(i+1, M), i] = t
        H_n[mod1(i-1, M), i] = t
    end
    return H_n
end

@doc """
Hamiltonian of isolated bar used for RGF.
"""
function isolated_bar(;M::Integer = 3, t::Real = 1.)
    H_n = spzeros(Float64, M^2, M^2)
    if M >= 2
        for m in 0:M-1, n in 0:M-1
            i = M*m + n + 1
            ip_x = M*mod(m+1, M) + n + 1
            im_x = M*mod(m-1, M) + n + 1
            ip_y = M*m + mod(n+1, M) + 1
            im_y = M*m + mod(n-1, M) + 1
            H_n[i, ip_x] = t
            H_n[i, im_x] = t
            H_n[i, ip_y] = t
            H_n[i, im_y] = t
            H_n[ip_x, i] = t
            H_n[im_x, i] = t
            H_n[ip_y, i] = t
            H_n[im_y, i] = t
        end

    end
    return H_n
end

function tmm_bar(; M::T = 3, E::F = 0., W::F = 1., t::F = 1., rng = Random.GLOBAL_RNG, N::T = 10^4, N_qr::T = 10) where {T <: Integer, F <: AbstractFloat}
    Msq = M^2
    if N%N_qr != 0
        error("N should be divisible by N_qr")
    end
    H_n = Matrix(isolated_bar(M = M, t = t))
    H_n = convert.(F, H_n)
    U = convert(Matrix{F}, I(2*Msq))
    U = U[:,1:Msq]
    ΣlnR = 0.
    IMsq = convert.(F, I(Msq))
    zerosq = zeros(F, Msq, Msq)

    EMsq = E*IMsq

    U_dummy = Matrix([zerosq IMsq; -IMsq zerosq])
    U_dummy_view = view(U_dummy, 1:Msq, 1:Msq)
    r_on = Array{F}(undef, Msq)
    for i in 1:N
        rand!(rng, r_on)
        r_on .-=  0.5
        rmul!(r_on, W)
        U_dummy_view .= (H_n .+ Diagonal(r_on) .- EMsq)
        U .= U_dummy*U
        if i%N_qr == 0
            U_new, R = qr!(U)
            Rs = sort(real.(R[diagind(R)]), rev = true)
            if i > 1000; display(Rs); end
            U .= @view Matrix(U_new)[:, 1:Msq]
            ΣlnR += log(abs(R[end, end]))
        end
    end
    return  N / ΣlnR #Full lypunov exponents
end
