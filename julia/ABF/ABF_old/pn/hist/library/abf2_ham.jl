# Construct/overwrite/modify full hamiltonian
# overwriting functions are useful for minimum memory access.
# But it must be used with existing variables or variable with zeros.

using Random
using LinearAlgebra
using StaticArrays
include("abf2_block.jl")


@doc """
Overwrite real hamiltonain created by block H0, H1, to the already existing nearest-neighbor hopping real hamiltonian
To create hamiltonian with this, H should be zeros square matrix
"""
function overwrite_ham!(H::AbstractArray{Float64, 2}, H0::AbstractArray{Float64, 2}, H1::AbstractArray{Float64, 2}, N::Int)
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
function create_ham(H0::AbstractArray{Float64, 2}, H1::AbstractArray{Float64, 2}, N::Int)
    H = zeros(Float64, 2N,2N)
    overwrite_ham!(H, H0, H1, N)
    return H
end

@doc """
Return uncorrelated onsite disorder to the input hamiltonian H.
W: disorder strength, H: matrix where onsite disorder is added, r_arr: array of random numbers
"""
function add_diag(H::AbstractArray{Float64, 2}, r_arr::AbstractVector{Float64})
    H_dis = copy(H)
    view(H_dis, diagind(H)) .+= r_arr
    return H_dis
end

@doc """
Add uncorrelated onsite disorder to hamiltonian H.
W: disorder strength, H: matrix where onsite disorder is added, r_arr: array of random numbers
"""
function add_diag!(H::AbstractArray{Float64, 2},r_arr::AbstractVector{Float64})
    broadcast!(+,view(H, diagind(H)), view(H, diagind(H)), r_arr)
    # for i = 1:size(H,1)
    #     H[i,i] += r_arr[i]
    # end
end

function sub_diag!(H::AbstractArray{Float64, 2},r_arr::AbstractVector{Float64})
    broadcast!(-,view(H, diagind(H)), view(H, diagind(H)), r_arr)
    # for i = 1:size(H,1)
    #     H[i,i] += r_arr[i]
    # end
end

@doc """
Make real ABF2 hamiltonian(E= -1, 1) with angles
unit of pi radian.
θ1, θ2: Local unitary transformation angles N: Number of unit cells.
"""
function create_ham_abf2(θ1::Real, θ2::Real, N::Int)
    U_1 = unitary(θ1)
    U_2 = unitary(θ2)
    H0, H1 = abf2_block(θ1, θ2);
    H = create_ham(H0, H1, N)
    return H
end

@doc """
Similar to create_H_ABF2, but overwrite for saving space.
"""
function overwrite_ham_abf2!(H::AbstractArray{Float64, 2}, θ1::Float64, θ2::Float64, N::Int)
    U_1 = unitary(θ1)
    U_2 = unitary(θ2)
    H0, H1 = abf2_block(θ1, θ2);
    overwrite_ham!(H, H0, H1, N)
end
