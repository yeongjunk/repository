using LinearAlgebra
using SparseArrays

export ham_fd, LUT, redef1, redef2, project!, project, U_fe, ham_fe

function ham_fd(ltc::Lattice2D, Ea::Real, Eb::Real)
    @assert ltc.U == 2
    return spdiagm(0 => repeat([Float64(Ea), Float64(Eb)], ltc.M*ltc.N))
end

function LUT(ltc::Lattice2D, θ::Real)
    num_sites = ltc.M*ltc.N*ltc.U
    cos_θ = cospi(θ)
    sin_θ = sinpi(θ)
    I = Int64[]; J = Int64[]; V = Float64[]
    for m in 1:ltc.M, n in 1:ltc.N
        push!(I, index(ltc, (m,n,1)))
        push!(J, index(ltc, (m,n,1)))
        push!(V, cos_θ )

        push!(I, index(ltc, (m,n,2)))
        push!(J, index(ltc, (m,n,1)))
        push!(V, -sin_θ )

        push!(I, index(ltc, (m,n,1)))
        push!(J, index(ltc, (m,n,2)))
        push!(V, sin_θ )

        push!(I, index(ltc, (m,n,2)))
        push!(J, index(ltc, (m,n,2)))
        push!(V, cos_θ )
    end
    return sparse(I, J, V, num_sites, num_sites)
end

function LUT(ltc::Lattice2D, θ::Real, ϕ1::Real, ϕ2::Real)
    num_sites = ltc.M*ltc.N*ltc.U
    cos_θ = cospi(θ)
    sin_θ = sinpi(θ)
    exp_phi1 = exp(im*ϕ1)
    exp_phi2 = exp(im*ϕ2)
    exp_phi1m = 1/exp_phi1
    exp_phi2m = 1/exp_phi2

    I = Int64[]; J = Int64[]; V = ComplexF64[]
    for m in 1:ltc.M, n in 1:ltc.N
        push!(I, index(ltc, (m,n,1)))
        push!(J, index(ltc, (m,n,1)))
        push!(V, exp_phi1*cos_θ)

        push!(I, index(ltc, (m,n,2)))
        push!(J, index(ltc, (m,n,1)))
        push!(V, -exp_phi2m*sin_θ)

        push!(I, index(ltc, (m,n,1)))
        push!(J, index(ltc, (m,n,2)))
        push!(V, exp_phi2*sin_θ)

        push!(I, index(ltc, (m,n,2)))
        push!(J, index(ltc, (m,n,2)))
        push!(V, exp_phi1m*cos_θ)
    end
    return sparse(I, J, V, num_sites, num_sites)
end

function redef1(ltc::Lattice2D)
    num_sites = ltc.M*ltc.N*ltc.U
    I = Int64[]; J = Int64[]; V = Float64[]
    for m in 1:ltc.M, n in 1:ltc.N
        push!(I, index(ltc, (m,n,1)))
        push!(J, index(ltc, (m,n,1)))
        push!(V, 1)

        push!(I, index(ltc, (m,n,2)))
        push!(J, index(ltc, (m,n+1,2)))
        push!(V, 1)
    end
    return sparse(I, J, V, num_sites, num_sites)
end

function redef2(ltc::Lattice2D)
    num_sites = ltc.M*ltc.N*ltc.U
    I = Int64[]; J = Int64[]; V = Float64[]
    for m in 1:ltc.M, n in 1:ltc.N
        push!(I, index(ltc, (m,n,1)))
        push!(J, index(ltc, (m,n,1)))
        push!(V, 1)

        push!(I, index(ltc, (m,n,2)))
        push!(J, index(ltc, (m+1,n,2)))
        push!(V, 1)
    end
    return sparse(I, J, V, num_sites, num_sites)
end

function U_fe(ltc::Lattice2D, θ::Real)
    U1 = LUT(ltc,θ)
    T1 = redef1(ltc)
    T2 = redef2(ltc)
    return U1*T2*U1*T1*U1
end

function U_fe(ltc::Lattice2D, θ::Real, ϕ1::Real, ϕ2::Real)
    U1 = LUT(ltc,θ, ϕ1, ϕ2)
    T1 = redef1(ltc)
    T2 = redef2(ltc)
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
