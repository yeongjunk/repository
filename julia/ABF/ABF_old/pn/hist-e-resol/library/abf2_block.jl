# Construct/Transform block matrices H0 , H1, U

using StaticArrays
using LinearAlgebra
@doc """
2 X 2 Real unitary matrix.
The unit of angle is pi(rad).
"""
function unitary(θ)
    c = cospi(θ)
    s = sinpi(θ)
    U = SA_F64[c s; -s c]
    return U
end

@doc """
transformation that results in unit cell redefinition (c_n -> a_n+1, d_n ->b_n)
to the H0, H1 blocks. only works when every H0 are same
"""
function T_p(H0::SMatrix{2, 2, Float64})
    T1 = @SMatrix [0. 0.; 0. 1.];
    T2 = @SMatrix [1. 0.; 0. 0.] # Matrices that changes the unit cell definition
    H0_new = T2*H0*T2 + T1*H0*T1
    H1_new = T2*H0*T1
    return H0_new, H1_new
end


@doc """
transformation that results in unit cell redefinition (c_n -> a_n-1, d_n ->b_n)
to the H0, H1 blocks (ignores H2, H3,... and so on)
"""
function T_n(H0::AbstractArray{Float64, 2}, H1::AbstractArray{Float64, 2})
    convert(SMatrix{2,2,Float64}, H0)
    convert(SMatrix{2,2,Float64}, H1) # conversion for typesafe operation

    T1 = @SMatrix [0. 0.;0. 1.];
    T2 = @SMatrix [1. 0.; 0. 0.] # Matrices that changes the unit cell definition
    H0_new = T2*H0*T2 + T1*H1'*T2 + T2*H1*T1 + T1*H0*T1
    H1_new = T2*H1*T2 + T1*H0*T2 + T1*H1*T1
    return H0_new, H1_new
end

@doc """
transformation that results in unit cell redefinition (c_n -> a_n-1, d_n ->b_n)
to the H0, H1 blocks (ignores H2, H3,... and so on). only works when all H0 and H1 are same
"""
function T_n(H0::SMatrix{2, 2, Float64}, H1::SMatrix{2, 2, Float64})
    T1 = @SMatrix [0. 0.;0. 1.];
    T2 = @SMatrix [1. 0.; 0. 0.] # Matrices that changes the unit cell definition
    H0_new = T2*H0*T2 + T1*H1'*T2 + T2*H1*T1 + T1*H0*T1
    H1_new = T2*H1*T2 + T1*H0*T2 + T1*H1*T1
    return H0_new, H1_new
end

@doc """
Construct the H0 and H1 of ABF2 real creutz ladder, E = -1,1 by applying U1 and U2
"""
function abf2_block(θ1::Float64, θ2::Float64)
    U_1 = unitary(θ1)
    U_2 = unitary(θ2)
    H0 = @SMatrix [1. 0.;0. -1.];

    H0 = U_1*H0*U_1'
    H0, H1 = T_p(H0)

    H0 = U_2*H0*U_2'
    H1 = U_2*H1*U_2'

    return H0, H1
end

@doc """
Inverse construction of the H0 and H1 of ABF2, E = -1,1 by applying U1 and U2
"""
function abf2_block_inv(θ1::Float64, θ2::Float64)
    U_1 = unitary(θ1)
    U_2 = unitary(θ2)

    H0 = @SMatrix [1. 0.;0. -1.]; H1 = @SMatrix [0. 0.;0. 0.]

    H1 = U_2'*H1*U_2
    H0 = U_2'*H0*U_2
    H0, H1 = T_n(H0, H1)
    H0 = U_1'*H0*U_1
    H1 = U_1'*H1*U_1


    return H0, H1
end
