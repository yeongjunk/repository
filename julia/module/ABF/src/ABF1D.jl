using LinearAlgebra
using SparseArrays

export ham_fd, LUT, redef1, project, U_fe, ham_fe

function ham_fd(ltc::Lattice1D, Ea::Real, Eb::Real)
    @assert ltc.U == 2
    return spdiagm(0 => repeat([Float64(Ea), Float64(Eb)], ltc.N))
end

function LUT(ltc::Lattice1D, θ::Real)
    num_sites = ltc.N*ltc.U
    cos_θ = cospi(θ)
    sin_θ = sinpi(θ)
    I = Int64[]; J = Int64[]; V = Float64[]
    for n in 1:ltc.N
        push!(I, index(ltc, (n, 1)))
        push!(J, index(ltc, (n, 1)))
        push!(V, cos_θ)

        push!(I, index(ltc, (n, 2)))
        push!(J, index(ltc, (n, 1)))
        push!(V, -sin_θ)

        push!(I, index(ltc, (n, 1)))
        push!(J, index(ltc, (n, 2)))
        push!(V, sin_θ)

        push!(I, index(ltc, (n, 2)))
        push!(J, index(ltc, (n, 2)))
        push!(V, cos_θ)
    end
    return sparse(I, J, V, num_sites, num_sites)
end

function redef1(ltc::Lattice1D)
    num_sites = ltc.N*ltc.U
    I = Int64[]; J = Int64[]; V = Float64[]
    for n in 1:ltc.N
        push!(I, index(ltc, (n, 1)))
        push!(J, index(ltc, (n, 1)))
        push!(V, 1)

        push!(I, index(ltc, (n, 2)))
        push!(J, index(ltc, (n + 1, 2)))
        push!(V, 1)
    end
    return sparse(I, J, V, num_sites, num_sites)
end

function U_fe(ltc::Lattice1D, θ::Real)
    U1 = LUT(ltc, θ)
    T1 = redef1(ltc)
    return U1*T1*U1
end

function ham_fe(ltc::Lattice1D, Ea::Real, Eb::Real, θ::Real)
    U = U_fe(ltc, θ)
    H_fd = ham_fd(ltc, Ea, Eb)
    H_fe = U*H_fd*U'
    return H_fe, U
end
