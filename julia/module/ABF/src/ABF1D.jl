using LinearAlgebra
using SparseArrays

export ham_fd, LUT, redef1, project, U_fe, ham_fe, ham_fe_obc, U_fe_obc, redef1_obc



"""
Creates fully detangled Hamiltonian of Float type F.
"""
function ham_fd(ltc::Lattice1D, Ea::F, Eb::F) where F
    @assert ltc.U == 2
    return spdiagm(0 => repeat([F(Ea), F(Eb)], ltc.N))
end


"""
Creates fully detangled Hamiltonian of Float type F, of size N and two band at Ea and Eb
"""
function ham_fd(N::Int64, Ea::F, Eb::F) where F
    ltc = Lattice1D(N, 2)
    return ham_fd(ltc, Ea, Eb) 
end

function LUT(ltc::Lattice1D, θ::F; pirad::Bool=true) where F 
    num_sites = ltc.N*ltc.U
    if pirad
        cos_θ = cospi(θ)
        sin_θ = sinpi(θ)
    else
        cos_θ = cos(θ)
        sin_θ = sin(θ)
    end
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


function LUT(N::Int64, θ::F; pirad::Bool=true) where F 
    ltc = Lattice1D(N, U)
    return LUT(ltc, θ, pirad=pirad)
end


function LUT(ltc::Lattice1D, θ::F, ϕ1::F, ϕ2::F) where F 
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


function LUT(N::Int64, θ::F, ϕ1::F, ϕ2::F) where F 
    ltc = Lattice1D(N, 2)
    return LUT(ltc, θ, ϕ1, ϕ2)
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
Unit cell redefinition transformation. default eltype is Float64. If you want more general type, specify keyword argument, e.g. vartype=BigFloat. 
"""
function redef1(N::Int64; vartype=Float64)
   ltc = Lattice1D(N, 2)
   return redef1(ltc, N, vartype=vartype) 
end

"""
Construct Full unitary that constructs FE ABF
"""
function U_fe(ltc::Lattice1D, θ::F; pirad=true) where F 
    U1 = LUT(ltc, θ, pirad=pirad)
    T1 = redef1(ltc, vartype = F)
    return U1*T1*U1
end


"""
Construct Full unitary that constructs FE ABF
"""
function U_fe(N::Int64, θ::F; pirad=true) where F 
    ltc = Lattice1D(N, 2)
    return U_fe(ltc, θ, pirad=pirad) 
end

"""
Real unitary parameters. This is the main FE construction
"""
function ham_fe(ltc::Lattice1D, Ea::F, Eb::F, θ::F; pirad=true) where F
    U = U_fe(ltc, θ, pirad=pirad)
    H_fd = ham_fd(ltc, Ea, Eb)
    H_fe = U*H_fd*U'
    return H_fe, U
end


"""
Real unitary parameters. This is the main FE construction
"""
function ham_fe(N::Int64, Ea::F, Eb::F, θ::F; pirad=true) where F
    ltc = Lattice1D(N, 2)
    return ham_fe(ltc, Ea, Eb, θ, pirad=pirad) 
end
