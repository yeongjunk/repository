using LinearAlgebra
using SparseArrays
using LinearMaps

export ham_fd, lut, redef1, U_fe, ham_fe

function ham_fd(ltc::Lattice1D, Ea::F, Eb::F) where F
    return LinearMap{F}((y, x) -> ham_fd!(y, x, Ea, Eb), 2ltc.N, ishermitian = true)
end

function lut(ltc::Lattice1D, θ::F; kwargs...) where F
    return LinearMap{F}((y, x) -> lut!(y, x, θ; kwargs...), (y, x) -> lut_adj!(y, x, θ; kwargs...), 2ltc.N)
end

"""
Unit cell redefinition transformation. default eltype is Float64. If you want more general type, specify keyword argument, e.g. vartype=BigFloat. 
"""
function redef1!(y::AbstractVector, x::AbstractVector, ltc::Lattice1D)
    for i in 1:2:ltc.N*ltc.U-2
        y[i]   = x[i]
        y[i+1] = x[i+3]
    end
    y[end-1] = x[end-1]
    y[end] = x[2]
end

"""
Unit cell redefinition transformation. default eltype is Float64. If you want more general type, specify keyword argument, e.g. vartype=BigFloat. 
"""
function redef1_adj!(y::AbstractVector, x::AbstractVector, ltc::Lattice1D)
    for i in 3:2:ltc.N*ltc.U
        y[i]   = x[i]
        y[i+1] = x[i-1]
    end
    y[1] = x[1]
    y[2] = x[end]
end

function redef1(ltc::Lattice1D)
    return LinearMap((y, x) -> redef1!(y, x, ltc), (y, x) -> redef1_adj!(y, x, ltc), 2ltc.N)
end

"""
Construct Full unitary that constructs FE ABF
"""
function U_fe(ltc::Lattice1D, θ::F; pirad::Bool=true) where F 
    U = lut(ltc, θ; pirad = pirad) 
    T = redef1(ltc) 
    return U*T*U 
end

"""
Real unitary parameters. This is the main FE construction
"""
function ham_fe(ltc::Lattice1D, Ea::F, Eb::F, θ::F; pirad=true) where F
    U = U_fe(ltc, θ, pirad=pirad)
    H_fd = ham_fd(ltc, Ea, Eb)
    H_fe = LinearMap{F}(U*H_fd*U', ishermitian=true)
    return H_fe, U
end
