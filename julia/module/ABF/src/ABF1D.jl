using LinearAlgebra
using SparseArrays

export ham_fd, LUT, redef1, project, U_fe, ham_fe

"""
Creates fully detangled Hamiltonian of Float type F.
"""
function ham_fd(ltc::Lattice1D, Ea::F, Eb::F) where F <: AbstractFloat
    @assert ltc.U == 2
    return spdiagm(0 => repeat([F(Ea), F(Eb)], ltc.N))
end

function LUT(ltc::Lattice1D, θ::F) where F <: AbstractFloat 
    num_sites = ltc.N*ltc.U
    cos_θ = cospi(θ)
    sin_θ = sinpi(θ)
    I = Int64[]; J = Int64[]; V = F[]
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

function LUT(ltc::Lattice1D, θ::F, ϕ1::F, ϕ2::F) where F <: AbstractFloat
    num_sites = ltc.N*ltc.U
    cos_θ = cospi(θ)
    sin_θ = sinpi(θ)
    exp_phi1 = exp(im*ϕ1)
    exp_phi2 = exp(im*ϕ2)
    exp_phi1m = 1/exp_phi1
    exp_phi2m = 1/exp_phi2

    I = Int64[]; J = Int64[]; V = Complex{F}[]
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


"""
Unit cell redefinition transformation. default eltype is Float64. If you want more general type, specify keyword argument, e.g. vartype=BigFloat. 
"""
function redef1(ltc::Lattice1D; vartype=Float64)
    num_sites = ltc.N*ltc.U
    I = Int64[]; J = Int64[]; V = vartype[]
    for n in 1:ltc.N
        push!(I, index(ltc, (n, 1)))
        push!(J, index(ltc, (n, 1)))
        push!(V, one(vartype))

        push!(I, index(ltc, (n, 2)))
        push!(J, index(ltc, (n + 1, 2)))
        push!(V, one(vartype))
    end
    return sparse(I, J, V, num_sites, num_sites)
end

"""
This is for test. Don't use it
"""
function redef1(ltc::Lattice1D, ϕ::F, vartype = F) where F <: AbstractFloat
    num_sites = ltc.N*ltc.U
    I = Int64[]; J = Int64[]; V = Complex{F}[]
    c = exp(im*ϕ)
    for n in 1:ltc.N
        push!(I, index(ltc, (n, 1)))
        push!(J, index(ltc, (n, 1)))
        push!(V, one(vartype))

        push!(I, index(ltc, (n, 2)))
        push!(J, index(ltc, (n + 1, 2)))
        push!(V, c)
    end
    return sparse(I, J, V, num_sites, num_sites)
end

"""
Construct Full unitary that constructs FE ABF
"""
function U_fe(ltc::Lattice1D, θ::F) where F <: AbstractFloat
    U1 = LUT(ltc, θ)
    T1 = redef1(ltc, vartype = F)
    return U1*T1*U1
end

"""
Real unitary parameters. This is the main FE construction
"""
function ham_fe(ltc::Lattice1D, Ea::F, Eb::F, θ::F) where F <: AbstractFloat
    U = U_fe(ltc, θ)
    H_fd = ham_fd(ltc, Ea, Eb)
    H_fe = U*H_fd*U'
    return H_fe, U
end

"""
More general complex LUT. You won't need it.
"""
function ham_fe(ltc::Lattice1D, Ea::F, Eb::F, θ::F, ϕ1::F, ϕ2::F) where F <: AbstractFloat
    U = U_fe(ltc, θ, ϕ1, ϕ2)
    H_fd = ham_fd(ltc, Ea, Eb)
    H_fe = U*H_fd*U'
    return H_fe, U
end
"""
Most general LUT parameters. You won't need it.
"""
function ham_fe(ltc::Lattice1D, Ea::F, Eb::F, θ::F, ϕ11::F, ϕ12::F, ϕ21::F, ϕ22::F, ϕ::F) where F <: AbstractFloat
    U1 = LUT(ltc, θ, ϕ11, ϕ12)
    T = redef1(ltc, ϕ)
    U2 = LUT(ltc, θ, ϕ21, ϕ22)
    U = U2*T*U1
    H_fd = ham_fd(ltc, Ea, Eb)
    H_fe = U*H_fd*U'
    return H_fe, U
end
