
using LinearAlgebra
using SparseArrays

export ham_fd, LUT, redef1, redef2, redef3, U_fe

function ham_fd(ltc::Lattice3D, Ea::F, Eb::F) where F <: AbstractFloat
    @assert ltc.U == 2
    return spdiagm(0 => repeat([F(Ea), F(Eb)], ltc.L*ltc.M*ltc.N))
end


function ham_fd(L::Integer, M::Integer, N::Integer, Ea::F, Eb::F) where F <: AbstractFloat
    ltc = Lattice3D(L, M, N, 2)
    return ham_fd(ltc, Ea, Eb) 
end

function LUT(ltc::Lattice3D, θ::F; pirad::Bool = true) where F <: AbstractFloat
    num_sites = ltc.L*ltc.M*ltc.N*ltc.U
    if pirad
        cos_θ = cospi(θ)
        sin_θ = sinpi(θ)
    else
        cos_θ = cos(θ)
        sin_θ = sin(θ)
    end

    I = Int64[]; J = Int64[]; V = F[]
    for l in 1:ltc.L, m in 1:ltc.M, n in 1:ltc.N
        push!(I, index(ltc, (l,m,n,1)))
        push!(J, index(ltc, (l,m,n,1)))
        push!(V, cos_θ )

        push!(I, index(ltc, (l,m,n,2)))
        push!(J, index(ltc, (l,m,n,1)))
        push!(V, -sin_θ )

        push!(I, index(ltc, (l,m,n,1)))
        push!(J, index(ltc, (l,m,n,2)))
        push!(V, sin_θ )

        push!(I, index(ltc, (l,m,n,2)))
        push!(J, index(ltc, (l,m,n,2)))
        push!(V, cos_θ )
    end
    return sparse(I, J, V, num_sites, num_sites)

end


function LUT(L::Integer, M::Integer, N::Integer, θ::F; pirad::Bool=true) where F <: AbstractFloat
    ltc = Lattice3D(L, M, N, 2)
    return LUT(ltc, θ, pirad=pirad)
end

function LUT(ltc::Lattice3D, θ::F, ϕ1::F, ϕ2::F; pirad=true) where F <: AbstractFloat
    num_sites = ltc.L*ltc.M*ltc.N*ltc.U
    if pirad
        cos_θ = cospi(θ)
        sin_θ = sinpi(θ)
    else
        cos_θ = cos(θ)
        sin_θ = sin(θ)
    end
    exp_phi1 = exp(im*ϕ1)
    exp_phi2 = exp(im*ϕ2)
    exp_phi1m = 1/exp_phi1
    exp_phi2m = 1/exp_phi2

    I = Int64[]; J = Int64[]; V = Complex{F}[]
    for l in 1:ltc.L, m in 1:ltc.M, n in 1:ltc.N
        push!(I, index(ltc, (l,m,n,1)))
        push!(J, index(ltc, (l,m,n,1)))
        push!(V, exp_phi1*cos_θ )

        push!(I, index(ltc, (l,m,n,2)))
        push!(J, index(ltc, (l,m,n,1)))
        push!(V, -exp_phi2m*sin_θ )

        push!(I, index(ltc, (l,m,n,1)))
        push!(J, index(ltc, (l,m,n,2)))
        push!(V, exp_phi2*sin_θ )

        push!(I, index(ltc, (l,m,n,2)))
        push!(J, index(ltc, (l,m,n,2)))
        push!(V, exp_phi1m*cos_θ )
    end
    return sparse(I, J, V, num_sites, num_sites)
end

function LUT(L::Integer, M::Integer, N::Integer, θ::F, ϕ1::F, ϕ2::F; pirad=true) where F
    ltc = Lattice3D(L, M, N, 2)
    return LUT(ltc, θ, ϕ1, ϕ2, pirad=pirad)
end

function redef1(ltc::Lattice3D; vartype = Float64)
    num_sites = ltc.L*ltc.M*ltc.N*ltc.U
    I = Int64[]; J = Int64[]; V = vartype[]
    for l in 1:ltc.L, m in 1:ltc.M, n in 1:ltc.N
        push!(I, index(ltc, (l,m,n,1)))
        push!(J, index(ltc, (l,m,n,1)))
        push!(V, 1)

        push!(I, index(ltc, (l,m,n,2)))
        push!(J, index(ltc, (l+1,m,n,2)))
        push!(V, 1)
    end
    return sparse(I, J, V, num_sites, num_sites)
end

function redef1(L::Integer, M::Integer, N::Integer; vartype = Float64)
    ltc = Lattice3D(L, M, N, 2)
    return redef1(ltc, vartype=vartype)
end

function redef2(ltc::Lattice3D; vartype = Float64)
    num_sites = ltc.L*ltc.M*ltc.N*ltc.U
    I = Int64[]; J = Int64[]; V = vartype[]
    for l in 1:ltc.L, m in 1:ltc.M, n in 1:ltc.N
        push!(I, index(ltc, (l,m,n,1)))
        push!(J, index(ltc, (l,m,n,1)))
        push!(V, 1)

        push!(I, index(ltc, (l,m,n,2)))
        push!(J, index(ltc, (l,m+1,n,2)))
        push!(V, 1)
    end
    return sparse(I, J, V, num_sites, num_sites)
end


function redef2(L::Integer, M::Integer, N::Integer; vartype = Float64)
    ltc = Lattice3D(L, M, N, 2)
    return redef2(ltc, vartype=vartype)
end

function redef3(ltc::Lattice3D; vartype = Float64)
    num_sites = ltc.L*ltc.M*ltc.N*ltc.U
    I = Int64[]; J = Int64[]; V = vartype[]
    for l in 1:ltc.L, m in 1:ltc.M, n in 1:ltc.N
        push!(I, index(ltc, (l,m,n,1)))
        push!(J, index(ltc, (l,m,n,1)))
        push!(V, 1)

        push!(I, index(ltc, (l,m,n,2)))
        push!(J, index(ltc, (l,m,n+1,2)))
        push!(V, 1)
    end
    return sparse(I, J, V, num_sites, num_sites)
end

function redef3(L::Integer, M::Integer, N::Integer; vartype = Float64)
    ltc = Lattice3D(L, M, N, 2)
    return redef3(ltc, vartype=vartype)
end

function U_fe(ltc::Lattice3D, θ::F; pirad=true) where F <: AbstractFloat
    U1 = LUT(ltc,θ, pirad=pirad)
    T1 = redef1(ltc, vartype = F)
    T2 = redef2(ltc, vartype = F)
    T3 = redef3(ltc, vartype = F)
    return U1*T3*U1*T2*U1*T1*U1
end


function U_fe(L::Integer, M::Integer, N::Integer, θ::F; pirad=true) where F <: AbstractFloat
    ltc = Lattice3D(L, M, N, 2)
    return U_fe(ltc, θ, pirad=pirad)
end

function ham_fe(ltc::Lattice3D, Ea::F, Eb::F, θ::F;pirad=true) where F <: AbstractFloat
    U = U_fe(ltc, θ, pirad=pirad)
    H_fd = ham_fd(ltc, Ea, Eb)
    H_fe = U*H_fd*U'
    return H_fe, U
end

function ham_fe(L::Integer, M::Integer, N::Integer, Ea::F, Eb::F, θ::F;pirad=true) where F <: AbstractFloat
    return ham_fe(Lattice3D(L, M, N, 2), Ea, Eb, θ, pirad=pirad) 
end

