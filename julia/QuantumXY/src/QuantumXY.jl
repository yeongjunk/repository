module QuantumXY

export σ_p, σ_m, σ_z, basis_vectors, ham_xy, ham_xy_kron,
    compute_magnetization, compute_correlation, basis_correlations, basis_magnetizations

using LinearAlgebra
using SparseArrays

const σ_p = sparse([0 1.; 0 0])
const σ_m = sparse([0 0; 1. 0])
const σ_z = sparse([1 0; 0 -1.])

@doc """
Full basis vectors encoded as up: 1, down: 0.
"""
function basis_vectors(; N = 3)
    basis = BitVector[]
    for i in 0:2^N-1
        x = BitVector(digits(UInt64(i), pad = N, base = 2))
        reverse!(x)
        pushfirst!(basis, x)
    end
    return basis
end

@doc """
Construct Hamiltonian of quantum XY spin model using kronecker products.
N: Number of spins
J: exchange interaction strength
h: external magnetic field strength
"""
function ham_xy_kron(;N = 3, J = 1., h = 1.)
    H = spzeros(Float64, 2^N, 2^N)
    for i = 0:N-2
        H .-= J*kron(I(2^i), σ_p, σ_m, I(2^(N-i-2)))
        H .-= h*kron(I(2^i), σ_z, I(2^(N-i-1)))
    end

    # Periodic boundary condition
    i = N-1
    H .-= (J*kron(σ_m, I(2^(N-2)), σ_p))'
    H .-= h*kron(I(2^i), σ_z)

    H .= Hermitian(H)
    dropzeros!(H)

    return H
end

@doc """
Construct Hamiltonian of quantum XY spin model using kronecker products.
N: Number of spins
J: exchange interaction strength
h: external magnetic field strength
"""
function ham_xy(; N = 3, J = 1., h = 1.)
    v = basis_vectors(N = N)
    down_up = BitVector([0; 1])
    row = Int64[]
    col = Int64[]
    val = Float64[]
    for i in 1:length(v)
        u = copy(v[i])
        onsite = 0.
        for ui in u
            ui ? onsite += h : onsite -= h
        end
        push!(row, i)
        push!(col, i)
        push!(val, -onsite)

        for j in 1:N-1
            if u[j:j+1] == down_up
                i_new = i - 2^(N-j-1)
                push!(row, i_new)
                push!(col, i)
                push!(val, -J)
            end
        end

        # Periodic boundary condition
        j = N
        if u[[j;1]] == down_up
            i_new = i + 2^(N-1) - 1
            push!(row, i)
            push!(col, i_new)
            push!(val, -J)
        end
    end
    H = dropzeros(sparse(Hermitian(sparse(row, col, val, 2^N, 2^N))))

    return H
end


function basis_magnetizations(; N = 3)
    basis = basis_vectors(N = N)
    m = zeros(Float64, 2^N)
    for (i, basis_i) in enumerate(basis)
        for s in basis_i
            s ? m[i] += 1 : m[i] -= 1
        end
    end
    return m/N
end

function basis_correlations(; N = 3, r = 1)
    basis = basis_vectors(N = N)
    corr = zeros(Float64, 2^N)
    m = zeros(Float64, 2^N)
    for (i, basis_i) in enumerate(basis)
        for s in basis_i
            s ? m[i] += 1 : m[i] -= 1
        end
    end

    for (i, basis_i) in enumerate(basis)
        for n in 1:N
            if !xor(basis_i[n], basis_i[mod1(n+r, N)])
                corr[i] += 1
            else
                basis = basis_vectors(N = N)
                corr[i] -= 1
            end
        end
    end
    return corr/N .- (m/N).^2
end

function compute_magnetization(v)
    N = convert(Int64, log2(length(v)))
    m_basis = basis_magnetizations(N = N)
    m = 0
    for i in 1:length(v)
        m += m_basis[i]*abs2(v[i])
    end
    return m
end

function compute_correlation(v, r)
    N = convert(Int64, log2(length(v)))
    corr_basis = basis_correlations(N = N, r = r)
    corr = 0
    for i in 1:length(v)
        corr += corr_basis[i]*abs2(v[i])
    end
    return corr
end

end # module
