using LinearAlgebra
using SparseArrays
using LinearMaps

export ham_fd, lut, redef1, U_fe, ham_fe, redef2

function ham_fd(ltc::Lattice2D, Ea::F, Eb::F) where F
    return LinearMap{F}((y, x) -> ham_fd!(y, x, Ea, Eb), 2ltc.N*ltc.M, ishermitian = true)
end

function lut(ltc::Lattice2D, θ::F; kwargs...) where F
    return LinearMap{F}((y, x) -> lut!(y, x, θ; kwargs...), (y, x) -> lut_adj!(y, x, θ; kwargs...), 2ltc.N*ltc.M)
end

"""
Unit cell redefinition transformation. default eltype is Float64. If you want more general type, specify keyword argument, e.g. vartype=BigFloat. 
"""
function redef1!(y::AbstractVector, x::AbstractVector, ltc::Lattice2D)
    for m in 1:ltc.M, n in 1:ltc.N 
        y[index(ltc, (m,n,1))] = x[index(ltc, (m,n,1))]
        y[index(ltc, (m,n,2))] = x[index(ltc, (m+1,n,2))]
    end
end

"""
Unit cell redefinition transformation. default eltype is Float64. If you want more general type, specify keyword argument, e.g. vartype=BigFloat. 
"""
function redef1_adj!(y::AbstractVector, x::AbstractVector, ltc::Lattice2D)
    for m in 1:ltc.M, n in 1:ltc.N 
        y[index(ltc, (m,n,1))] = x[index(ltc, (m,n,1))]
        y[index(ltc, (m,n,2))] = x[index(ltc, (m-1,n,2))]
    end
end

function redef1(ltc::Lattice2D)
    return LinearMap((y, x) -> redef1!(y, x, ltc), (y, x) -> redef1_adj!(y, x, ltc), 2ltc.M*ltc.N)
end


"""
Unit cell redefinition transformation. default eltype is Float64. If you want more general type, specify keyword argument, e.g. vartype=BigFloat. 
"""
function redef2!(y::AbstractVector, x::AbstractVector, ltc::Lattice2D)
    for m in 1:ltc.M, n in 1:ltc.N 
        y[index(ltc, (m,n,1))] = x[index(ltc, (m,n,1))]
        y[index(ltc, (m,n,2))] = x[index(ltc, (m,n+1,2))]
    end
end

"""
Unit cell redefinition transformation. default eltype is Float64. If you want more general type, specify keyword argument, e.g. vartype=BigFloat. 
"""
function redef2_adj!(y::AbstractVector, x::AbstractVector, ltc::Lattice2D)
    for m in 1:ltc.M, n in 1:ltc.N 
        y[index(ltc, (m,n,1))] = x[index(ltc, (m,n,1))]
        y[index(ltc, (m,n,2))] = x[index(ltc, (m,n-1,2))]
    end
end

function redef2(ltc::Lattice2D)
    return LinearMap((y, x) -> redef2!(y, x, ltc), (y, x) -> redef2_adj!(y, x, ltc), 2ltc.M*ltc.N)
end
"""
Construct Full unitary that constructs FE ABF
"""
function U_fe(ltc::Lattice2D, θ::F; pirad::Bool=true) where F 
    U = lut(ltc, θ; pirad = pirad) 
    Tx = redef1(ltc) 
    Ty = redef2(ltc) 
    return U*Tx*U*Ty*U 
end

"""
Real unitary parameters. This is the main FE construction
"""
function ham_fe(ltc::Lattice2D, Ea::F, Eb::F, θ::F; pirad=true) where F
    U = U_fe(ltc, θ, pirad=pirad)
    H_fd = ham_fd(ltc, Ea, Eb)
    H_fe = LinearMap{F}(U*H_fd*U', ishermitian=true)
    return H_fe, U
end
