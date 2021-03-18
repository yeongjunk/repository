# Matrix transformations that are capable of detangling/entangling operations
# FD(& effective mdoel) <-> SD <-> FE in nu = 2, d = 1 ABF lattice
# the blockdiagonal function returns the sparse matrix.
# This reduces the memory allocation dramatically

# toto: 1. generalization to higher dimensions.
# 2. matrix index <-> lattice position(braket)
# 3. Optimization

using LinearAlgebra
using SparseArrays

@doc """
Construct (mN X mN) sparse block diagonal matrix from (m X m) block A.
"""
function blockdiagonal(A_block, N)
    m = size(A_block,1) #estimate the size of block

    # Initialize a block matrix and fill with zeros
    rowidx = repeat(1:m*N, inner = [m])
    colidx = repeat(1:m:m*N, inner = [m*m]) + repeat(0:m-1, outer = [m*N])
    vals = repeat(reshape(transpose(A_block), m*m), N)

    A = sparse(rowidx, colidx, vals)
    return A
end

@doc """
Create 2 X 2 real LUT block. Angle is in unit of pi radian.
"""
function LUT(θ)
    c = cospi(θ)
    s = sinpi(θ)
    U = Matrix{ComplexF64}([c s; -s c])

    return U
end

@doc """
Construct 2N X 2N real LUT block diagonal matrix. Angle is in unit of pi radian.
"""
function LUT(θ, N)
    U_block = LUT(θ)
    U = blockdiagonal(U_block, N)
    return U
end


@doc """
Create 2N X 2N unit cell redefining sparse matrix (a_n, b_n) -> (a_(n+1) b_(n))
"""
function uc_redef(N)
    rowidx = collect(1:2N)
    colidx = collect(1:2N)
    colidx[1:2:2N-1] = circshift(colidx[1:2:2N-1],-1)
    vals = ones(ComplexF64, 2N)
    return sparse(rowidx, colidx, vals)
end


@doc """
Create 2N X 2N blcok diagonal detangled hamiltonain of energy gap Δ_FB.
Lower site(b_n), which corresponds to even indices, has the lower band energy.
"""
function ham_FD(Δ_FB, N)
    H_block = Matrix{ComplexF64}([1 0; 0 0])
    H_FD = blockdiagonal(Δ_FB*H_block, N)
    return H_FD
end


@doc """
Remove every odd elements from 2N x 2N hamiltonian, i.e. project to N x N hamiltonian.
"""
function projection(H_FD)
    N = size(H_FD, 1)÷2
    I, J, V = findnz(H_FD)

    idx = findall(iseven.(I) .& iseven.(J)) # Find all even
    return sparse(I[idx].÷2, J[idx].÷2, V[idx], N, N)
end
