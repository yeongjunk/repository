# ν = 2 local change of basis operations on wavefunctions


function splitpsi(psi::AbstractVector{Float64})
    n = length(psi)
    a = Array{Float64}(undef, div(n,2))
    b = Array{Float64}(undef, div(n,2))

    for i = 1:n÷2
        a[i] = psi[2*i-1]
        b[i] = psi[2*i]
    end

    return a, b
end

@doc """
unitcell redefinition transformation (a_n, b_n) -> (a_n+1,b_n)
"""
function uc_redef_p(a::AbstractVector{Float64},b::AbstractVector{Float64})
    c = circshift(a, -1)
    return c, b
end

@doc """
unitcell redefinition transformation (a_n, b_n) -> (a_n-1,b_n)
"""
function uc_redef_n(a::AbstractVector{Float64},b::AbstractVector{Float64})
    c = circshift(a, 1)
    return c, b
end

@doc """
local unitary transformation
"""
function unitary(a::AbstractVector{Float64},b::AbstractVector{Float64}, U::AbstractArray{Float64, 2})
    c = U[1,1]*a .+ U[1,2]*b
    d = U[2,1]*a .+ U[2,2]*b

    return c, d
end
