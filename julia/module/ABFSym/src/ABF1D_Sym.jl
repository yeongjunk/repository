using LinearAlgebra
using SparseArrays
using SymPy

export ham_fd, LUT, redef1, project, U_fe, ham_fe

function ham_fd(ltc::Lattice1D, Ea::F, Eb::F) where F
    @assert ltc.U == 4
    return spdiagm(0 => repeat([F(Ea), F(Ea), F(Eb), F(Eb)], ltc.N))
end

function LUT(ltc::Lattice1D, θ::F) where F
    num_sites = ltc.N*ltc.U
    cos_θ = cospi(θ)
    sin_θ = sinpi(θ)
    I = Int64[]; J = Int64[]; V = F[]
    for n in 1:ltc.N
        #Spin up 
        push!(I, index(ltc, (n, 1)))
        push!(J, index(ltc, (n, 1)))
        push!(V, cos_θ)

        #Spin down 
        push!(I, index(ltc, (n, 2)))
        push!(J, index(ltc, (n, 2)))
        push!(V, cos_θ)

        #Spin up 
        push!(I, index(ltc, (n, 3)))
        push!(J, index(ltc, (n, 3)))
        push!(V, cos_θ)

        #Spin down
        push!(I, index(ltc, (n, 4)))
        push!(J, index(ltc, (n, 4)))
        push!(V, cos_θ)

        #Spin up
        push!(I, index(ltc, (n, 1)))
        push!(J, index(ltc, (n, 3)))
        push!(V, sin_θ)

        #Spin down
        push!(I, index(ltc, (n, 2)))
        push!(J, index(ltc, (n, 4)))
        push!(V, sin_θ)

        #Spin up
        push!(I, index(ltc, (n, 3)))
        push!(J, index(ltc, (n, 1)))
        push!(V, -sin_θ)

        #Spin down
        push!(I, index(ltc, (n, 4)))
        push!(J, index(ltc, (n, 2)))
        push!(V, -sin_θ)

    end
    return sparse(I, J, V, num_sites, num_sites)
end


function LUT(ltc::Lattice1D, θ::Sym)
    num_sites = ltc.N*ltc.U
    cos_θ = cos(θ)
    sin_θ = sin(θ)
    I = Int64[]; J = Int64[]; V = Sym[]
    for n in 1:ltc.N
        #Spin up 
        push!(I, index(ltc, (n, 1)))
        push!(J, index(ltc, (n, 1)))
        push!(V, cos_θ)

        #Spin down 
        push!(I, index(ltc, (n, 2)))
        push!(J, index(ltc, (n, 2)))
        push!(V, cos_θ)

        #Spin up 
        push!(I, index(ltc, (n, 3)))
        push!(J, index(ltc, (n, 3)))
        push!(V, cos_θ)

        #Spin down
        push!(I, index(ltc, (n, 4)))
        push!(J, index(ltc, (n, 4)))
        push!(V, cos_θ)

        #Spin up
        push!(I, index(ltc, (n, 1)))
        push!(J, index(ltc, (n, 3)))
        push!(V, sin_θ)

        #Spin down
        push!(I, index(ltc, (n, 2)))
        push!(J, index(ltc, (n, 4)))
        push!(V, sin_θ)

        #Spin up
        push!(I, index(ltc, (n, 3)))
        push!(J, index(ltc, (n, 1)))
        push!(V, -sin_θ)

        #Spin down
        push!(I, index(ltc, (n, 4)))
        push!(J, index(ltc, (n, 2)))
        push!(V, -sin_θ)

    end
    return sparse(I, J, V, num_sites, num_sites)
end


function redef1(ltc::Lattice1D; elt = Float64)
    num_sites = ltc.N*ltc.U
    I = Int64[]; J = Int64[]; V = elt[]
    for n in 1:ltc.N
        #spin downs
        push!(I, index(ltc, (n, 1)))
        push!(J, index(ltc, (n, 1)))
        push!(V, 1)

        push!(I, index(ltc, (n, 3)))
        push!(J, index(ltc, (n + 1, 3)))
        push!(V, 1)

        #spin ups
        push!(I, index(ltc, (n, 2)))
        push!(J, index(ltc, (n, 2)))
        push!(V, 1)

        push!(I, index(ltc, (n, 4)))
        push!(J, index(ltc, (n + 1, 4)))
        push!(V, 1)
    end
    return sparse(I, J, V, num_sites, num_sites)
end

function U_fe(ltc::Lattice1D, θ::F) where F
    U1 = LUT(ltc, θ)
    T1 = redef1(ltc; elt = eltype(θ))
    return U1*T1*U1
end

function ham_fe(ltc::Lattice1D, Ea, Eb, θ)
    U = U_fe(ltc, θ)
    H_fd = ham_fd(ltc, Ea, Eb)
    H_fe = U*H_fd*U'
    return H_fe, U
end
