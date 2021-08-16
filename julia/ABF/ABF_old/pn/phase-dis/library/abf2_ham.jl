# Construct/overwrite/modify full hamiltonian
# overwriting functions are useful for minimum memory access.
# But it must be used with existing variables or variable with zeros.

using LinearAlgebra
include("abf2_block.jl")


@doc """
Overwrite real hamiltonain created by block H0, H1, to the already existing nearest-neighbor hopping real hamiltonian
To create hamiltonian with this, H should be zeros square matrix
"""
function overwrite_ham!(H, H0, H1, N)
    ν = 2
    M = ν*N

    if size(H) != (M, M) && (size(H0) && size(H1)) != (ν,ν)
        error("size error")
    end

    A = hcat(H1',H0,H1)

    ref = range(1, M, step = ν)
    for i in 2:length(ref)-1
        H[ref[i]:(ref[i+1]-1), ref[i-1]:(ref[i+1]+ν-1)] .= A
    end

    # seperately considered part to satisfy the periodic boundary condition
    H[ref[1]:ref[2]-1, ref[1]:ref[2]+ν-1] .= hcat(H0,H1)
    H[ref[1]:ref[2]-1, mod1(ref[1]-ν,M):mod1(ref[1]-1,M)] .= H1'

    H[ref[end]:end, ref[end-1]:ref[end]+ν-1] .= hcat(H1',H0)
    H[ref[end]:end, mod1(ref[end]+ν, M):mod1(ref[end]+2ν-1, M)] .= H1
end


@doc """
create a real hamiltonian created by blocks H0 and H1 with unit cell number N
"""
function create_ham(H0, H1, N)
    H = zeros(Float64, 2N,2N)
    overwrite_ham!(H, H0, H1, N)
    return H
end

@doc """
Add uncorrelated onsite disorder to hamiltonian H.
W: disorder strength, H: matrix where onsite disorder is added, r_on: array of random numbers
"""
function add_diag!(H, r_on)
    broadcast!(+, view(H, diagind(H)), view(H, diagind(H)), r_on)
    # for i = 1:size(H,1)
    #     H[i,i] += r_arr[i]
    # end
end

function add_diag(H,r_on)
    H_dis = copy(H)
    add_diag!(H_dis,r_on)
    return H_dis
end


function add_hop_dis!(H, r_hop)
    N = size(H,1) ÷ 2
    for j = 1:2N
        for i in 1:3
            H[j, mod1(j+i, 2N)] = H[j, mod1(j+i, 2N)]*exp(im*π*(r_hop[j,i]))
        end
    end

    H[1:2, end-1:end] = H[end-1:end,1:2]' #Periodic boundary Condition
    H = UpperTriangular(H) + UnitUpperTriangular(H)' - convert(Array{ComplexF64,2}, I(2N))
end

@doc """
Make real ABF2 hamiltonian(E= -1, 1) with angles
unit of pi radian.
θ1, θ2: Local unitary transformation angles N: Number of unit cells.
"""
function create_ham_abf2(θ1, θ2, N)
    U_1 = unitary(θ1)
    U_2 = unitary(θ2)
    H0, H1 = abf2_block(θ1, θ2);
    H = create_ham(H0, H1, N)
    return H
end

@doc """
Similar to create_H_ABF2, but overwrite for saving space.
"""
function overwrite_ham_abf2!(H, θ1, θ2, N)
    U_1 = unitary(θ1)
    U_2 = unitary(θ2)
    H0, H1 = abf2_block(θ1, θ2);
    overwrite_ham!(H, H0, H1, N)
end
