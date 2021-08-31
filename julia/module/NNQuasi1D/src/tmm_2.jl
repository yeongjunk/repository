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


@Threads.threads for i in 1:10
    println(i)
end

function tmm_bar(; M::T = 3, E::F = 0., W::F = 1., t::F = 1., rng = Random.GLOBAL_RNG, N::T = 10^4, q::T = 5, r::T = 100) where {T <: Integer, F <: AbstractFloat}
    Msq = M^2
    if N%q != 0
        error("N should be divisible by q")
    end
    if N%(q*r) != 0
        error("N should be divisible by q*r")
    end

    H_n = Matrix(isolated_bar(M = M, t = t))
    H_n = convert.(F, H_n)
    U = convert(Matrix{F}, I(2*Msq))
    U = U
    IMsq = I(Msq)
    EMsq = E*IMsq
    zerosq = zeros(F, Msq, Msq)

    U_dummy = Matrix{F}([zerosq IMsq; -IMsq zerosq])
    U_dummy_view = view(U_dummy, 1:Msq, 1:Msq)
    r_on = Vector{F}(undef, Msq)
    p = q*r
    s = N÷p

    D = 0.; Dsq = 0.; γ = 0.; γsq = 0.
    x = 0
    for i in 1:N
        rand!(rng, r_on)
        r_on .-=  0.5
        rmul!(r_on, W)
        U_dummy_view .= (H_n .+ Diagonal(r_on) .- EMsq)

        U .= U_dummy*U
        if i%q == 0
            U_new, R = qr!(U)
            if !issorted(abs.(R[diagind(R)]), rev = true)
                x += 1
            end
            display(R[diagind(R)])
            @views U .= Matrix{F}(U_new)
            log_R_i = log(minimum(abs.(R[diagind(R)])))
            D += log_R_i
            Dsq += log_R_i^2
        end

        if i%p == 0
            γ += D/p
            γsq += Dsq/p
            D = 0.
            Dsq = 0. #Initialize D
        end
    end
    println(x)
    γ /= s
    γsq /= s
    γ_stde = sqrt((γsq - γ^2)/(s - 1)/s)
    return  1/γ , γ_stde #Full lypunov exponents
end

function iter_mean_std(x)
    M = x[1]
    S = 0.
    @inbounds for k in 2:length(x)
        M_temp = M
        M += (x[k] - M)/k
        S += (x[k] - M_temp)*(x[k] - M)
    end
    return M, sqrt(S/(length(x)-1))
end
