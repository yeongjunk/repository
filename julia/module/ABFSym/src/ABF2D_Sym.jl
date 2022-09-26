using LinearAlgebra
using SparseArrays
using Lattices
export ham_fd, LUT, redef1, project, U_fe, ham_fe

function compute_num_sites(ltc::Lattice2D)
    return ltc.M*ltc.N*ltc.U
end
function ham_fd(ltc::Lattice2D, Ea::Real, Eb::Real)
    @assert ltc.U == 4
    return spdiagm(0 => repeat([Float64(Ea), Float64(Ea), Float64(Eb), Float64(Eb)], ltc.M*ltc.N))
end

function LUT(ltc::Lattice2D, θ::Real)
    num_sites = compute_num_sites(ltc)
    cos_θ = cospi(θ)
    sin_θ = sinpi(θ)
    I = Int64[]; J = Int64[]; V = Float64[]
    a↓ = 1
    a↑ = 2
    b↓ = 3
    b↑ = 4
    for m in 1:ltc.M, n in 1:ltc.N
        #Spin downs
        push!(I, index(ltc, (m, n, a↓)))
        push!(J, index(ltc, (m, n, a↓)))
        push!(V, cos_θ)

        push!(I, index(ltc, (m, n, b↓)))
        push!(J, index(ltc, (m, n, a↓)))
        push!(V, -sin_θ)

        push!(I, index(ltc, (m, n, a↓)))
        push!(J, index(ltc, (m, n, b↓)))
        push!(V, sin_θ)

        push!(I, index(ltc, (m, n, b↓)))
        push!(J, index(ltc, (m, n, b↓)))
        push!(V, cos_θ)
        #Spin ups
        push!(I, index(ltc, (m, n, a↑)))
        push!(J, index(ltc, (m, n, a↑)))
        push!(V, cos_θ)

        push!(I, index(ltc, (m, n, b↑)))
        push!(J, index(ltc, (m, n, a↑)))
        push!(V, -sin_θ)

        push!(I, index(ltc, (m, n, a↑)))
        push!(J, index(ltc, (m, n, b↑)))
        push!(V, sin_θ)

        push!(I, index(ltc, (m, n, b↑)))
        push!(J, index(ltc, (m, n, b↑)))
        push!(V, cos_θ)
    end
    return sparse(I, J, V, num_sites, num_sites)
end

function redef1(ltc::Lattice2D, ϕ::Real)
    num_sites = compute_num_sites(ltc)
    I = Int64[]; J = Int64[]; V = ComplexF64[]
    a↓ = 1
    a↑ = 2
    b↓ = 3
    b↑ = 4
    for m in 1:ltc.M, n in 1:ltc.N
        #spin downs
        push!(I, index(ltc, (m, n, a↓)))
        push!(J, index(ltc, (m, n, a↓)))
        push!(V, 1)

        push!(I, index(ltc, (m, n, b↓)))
        push!(J, index(ltc, (m, n + 1, b↓)))
        push!(V, exp(im*ϕ))

        #spin ups
        push!(I, index(ltc, (m, n, a↑)))
        push!(J, index(ltc, (m, n, a↑)))
        push!(V, 1)

        push!(I, index(ltc, (m, n, b↑)))
        push!(J, index(ltc, (m, n + 1, b↑)))
        push!(V, exp(im*ϕ))
    end
    return sparse(I, J, V, num_sites, num_sites)
end

function redef1(ltc::Lattice2D)
    return redef1(ltc, 0.)
end



function redef2(ltc::Lattice2D, ϕ::Real)
    num_sites = compute_num_sites(ltc)
    I = Int64[]; J = Int64[]; V = ComplexF64[]
    a↓ = 1
    a↑ = 2
    b↓ = 3
    b↑ = 4
    for m in 1:ltc.M, n in 1:ltc.N
        #spin downs
        push!(I, index(ltc, (m, n, a↓)))
        push!(J, index(ltc, (m, n, a↓)))
        push!(V, 1)

        push!(I, index(ltc, (m, n, b↓)))
        push!(J, index(ltc, (m + 1, n, b↓)))
        push!(V, exp(im*ϕ))

        #spin ups
        push!(I, index(ltc, (m, n, a↑)))
        push!(J, index(ltc, (m, n, a↑)))
        push!(V, 1)

        push!(I, index(ltc, (m, n, b↑)))
        push!(J, index(ltc, (m + 1, n, b↑)))
        push!(V, exp(im*ϕ))
    end
    return sparse(I, J, V, num_sites, num_sites)
end

function redef2(ltc::Lattice2D)
    return redef2(ltc, 0.)
end

function U_fe(ltc::Lattice2D, θ::Real)
    U1 = LUT(ltc, θ)
    T1 = redef1(ltc)
    T2 = redef2(ltc)
    return U1*T2*U1*T1*U1
end

function U_fe(ltc::Lattice2D, θ::Real, ϕ1::Real, ϕ2::Real)
    U1 = LUT(ltc, θ)
    T1 = redef1(ltc, ϕ1)
    T2 = redef2(ltc, ϕ2)
    return U1*T2*U1*T1*U1
end

function ham_fe(ltc::Lattice2D, Ea::Real, Eb::Real, θ::Real)
    U = U_fe(ltc, θ)
    H_fd = ham_fd(ltc, Ea, Eb)
    H_fe = U*H_fd*U'
    return H_fe, U
end

function ham_fe(ltc::Lattice2D, Ea::Real, Eb::Real, θ::Real, ϕ1::Real, ϕ2::Real)
    U = U_fe(ltc, θ, ϕ1, ϕ2)
    H_fd = ham_fd(ltc, Ea, Eb)
    H_fe = U*H_fd*U'
    return H_fe, U
end
