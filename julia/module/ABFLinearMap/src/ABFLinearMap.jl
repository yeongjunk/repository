module ABFLinearMap
    import Lattices.Lattice1D, Lattices.Lattice2D, Lattices.Lattice3D,
        Lattices.index, Lattices.site
    export ham_fd!, lut!, lut_adj!
    using LinearAlgebra
    using SparseArrays

    include("ABF1D.jl")
    include("ABF2D.jl")
#     include("ABF3D.jl")

"""
Creates fully detangled Hamiltonian of Float type F.
"""
function ham_fd!(y::AbstractVector, x::AbstractVector, Ea::F, Eb::F) where F
    for i in 1:2:length(x)
        y[i]   = Ea*x[i]
        y[i+1] = Eb*x[i+1] 
    end
    return y 
end

function lut!(y::AbstractVector, x::AbstractVector, θ::F; pirad::Bool=true) where F 
    if pirad
        cos_θ = cospi(θ)
        sin_θ = sinpi(θ)
    else
        cos_θ = cos(θ)
        sin_θ = sin(θ)
    end
    for i in 1:2:length(x)
        y[i]   =  cos_θ*x[i] + sin_θ*x[i+1]
        y[i+1] = -sin_θ*x[i] + cos_θ*x[i+1]
    end
    return y
end

function lut_adj!(y::AbstractVector, x::AbstractVector, θ::F; pirad::Bool=true) where F 
    if pirad
        cos_θ = cospi(θ)
        sin_θ = sinpi(θ)
    else
        cos_θ = cos(θ)
        sin_θ = sin(θ)
    end
    for i in 1:2:length(x)
        y[i]   =  cos_θ*x[i] - sin_θ*x[i+1]
        y[i+1] =  sin_θ*x[i] + cos_θ*x[i+1]
    end
    return y
end

    function project!(H_p::AbstractArray,H::AbstractArray)
        H_p .= H[2:2:end, 2:2:end]
    end

    function project(H::AbstractArray)
        H_p = H[2:2:end, 2:2:end]
        return H_p
    end

    export project!, project, Lattices
end
