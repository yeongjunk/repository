# Construct/Transform block matrices H0 , H1, U
using LinearAlgebra
using BlockArrays

@doc """
3 X 3 Real unitary matrix.
The unit of angle is pi(rad).
"""
function unitary_block(θ)::AbstractArray{ComplexF64}
    c = cospi(θ)
    s = sinpi(θ)
    U_bc = Array{ComplexF64}([1 0 0; 0 c -s; 0 s c])
    U_ab = Array{ComplexF64}([c -s 0; s c 0; 0 0 1])
    U = U_bc*U_ab
    return U
end

@doc """
A full mN X mN block diagonal matrix.
A is the input matrix.
"""
function blockdiagonal(A, N)
    m = size(A,1)

    # Initialize a block matrix and fill with zeros
    A_full = BlockArray{ComplexF64}(undef, repeat([m],N), repeat([m],N))
    fill!(A_full,0)
    for i in 1:N
        setblock!(A_full, A, i,i)
    end
    return A_full
end

@doc """
Redefine unit cell as (a'_n, b'_n, c'_n)=(a_n+1, b_n, c_n).
N is the number of unit cell.
"""
function uc_redef(N)
    T = Matrix{ComplexF64}(I, 3N, 3N)
    temp = T[1,:]
    for i in range(1,(3N-5), step = 3)
        T[i,:] = T[i+3,:]
    end
    T[3N-2,:] = temp
    return T
end
