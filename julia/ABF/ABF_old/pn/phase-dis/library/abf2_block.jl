# Construct/Transform block matrices H0 , H1, U

@doc """
2 X 2 Real unitary matrix.
The unit of angle is pi(rad).
"""
function unitary(θ)
    c = cospi(θ)
    s = sinpi(θ)
    U = [c s; -s c]
    return U
end

@doc """
transformation that results in unit cell redefinition (c_n -> a_n+1, d_n ->b_n)
to the H0, H1 blocks. only works when every H0 are same
"""
function T_p(H0)
    T1 = [0. 0.; 0. 1.];
    T2 = [1. 0.; 0. 0.]; # Matrices that changes the unit cell definition
    H0_new = T2*H0*T2 + T1*H0*T1
    H1_new = T2*H0*T1
    return H0_new, H1_new
end

@doc """
Construct the H0 and H1 of ABF2 real creutz ladder, E = -1,1 by applying U1 and U2
"""
function abf2_block(θ1, θ2)
    U_1 = unitary(θ1)
    U_2 = unitary(θ2)
    H0 = [1. 0.;0. -1.];

    H0 = U_1*H0*U_1'
    H0, H1 = T_p(H0)

    H0 = U_2*H0*U_2'
    H1 = U_2*H1*U_2'

    return H0, H1
end
