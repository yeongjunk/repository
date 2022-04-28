using LinearAlgebra
using SparseArrays

export ham_fd, LUT, redef1, project, U_fe, ham_fe

function ham_fd(ltc::Lattice1D, Ea::Real, Eb::Real)
    @assert ltc.U == 2
    return spdiagm(0 => repeat([Ea, Eb], ltc.N))
end

function LUT(ltc::Lattice1D, θ::Real)
    num_sites = ltc.N*ltc.U
    cos_θ = cospi(θ)
    sin_θ = sinpi(θ)
    I = Int64[]; J = Int64[]; V = eltype(θ)[]
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

function LUT(ltc::Lattice1D, θ::Real, ϕ1::Real, ϕ2::Real)
    num_sites = ltc.N*ltc.U
    cos_θ = cospi(θ)
    sin_θ = sinpi(θ)
    exp_phi1 = exp(im*ϕ1)
    exp_phi2 = exp(im*ϕ2)
    exp_phi1m = 1/exp_phi1
    exp_phi2m = 1/exp_phi2

    I = Int64[]; J = Int64[]; V = Complex{eltype(θ)}[]
    for n in 1:ltc.N
        push!(I, index(ltc, (n, 1)))
        push!(J, index(ltc, (n, 1)))
        push!(V, exp_phi1*cos_θ)

        push!(I, index(ltc, (n, 2)))
        push!(J, index(ltc, (n, 1)))
        push!(V, -exp_phi2*sin_θ)

        push!(I, index(ltc, (n, 1)))
        push!(J, index(ltc, (n, 2)))
        push!(V, exp_phi2m*sin_θ)

        push!(I, index(ltc, (n, 2)))
        push!(J, index(ltc, (n, 2)))
        push!(V, exp_phi1m*cos_θ)
    end
    return sparse(I, J, V, num_sites, num_sites)
end

function redef1(ltc::Lattice1D; tp = Float64)
    num_sites = ltc.N*ltc.U
    I = Int64[]; J = Int64[]; V = tp[]
    for n in 1:ltc.N
        push!(I, index(ltc, (n, 1)))
        push!(J, index(ltc, (n, 1)))
        push!(V, one(tp))

        push!(I, index(ltc, (n, 2)))
        push!(J, index(ltc, (n + 1, 2)))
        push!(V, one(tp))
    end
    return sparse(I, J, V, num_sites, num_sites)
end

function redef1(ltc::Lattice1D, ϕ::Real)
    num_sites = ltc.N*ltc.U
    I = Int64[]; J = Int64[]; V = eltype(ϕ)[]
    c = exp(im*ϕ)
    for n in 1:ltc.N
        push!(I, index(ltc, (n, 1)))
        push!(J, index(ltc, (n, 1)))
        push!(V, 1)

        push!(I, index(ltc, (n, 2)))
        push!(J, index(ltc, (n + 1, 2)))
        push!(V, c)
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

function ham_fe(ltc::Lattice1D, Ea::Real, Eb::Real, θ::Real, ϕ1::Real, ϕ2::Real)
    U = U_fe(ltc, θ, ϕ1, ϕ2)
    H_fd = ham_fd(ltc, Ea, Eb)
    H_fe = U*H_fd*U'
    return H_fe, U
end

function ham_fe(ltc::Lattice1D, Ea::Real, Eb::Real, θ::Real, ϕ11::Real, ϕ12::Real, ϕ21::Real, ϕ22::Real, ϕ::Real)
    U1 = LUT(ltc, θ, ϕ11, ϕ12)
    T = redef1(ltc, ϕ)
    U2 = LUT(ltc, θ, ϕ21, ϕ22)
    U = U2*T*U1
    H_fd = ham_fd(ltc, Ea, Eb)
    H_fe = U*H_fd*U'
    return H_fe, U
end
