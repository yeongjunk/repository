# Box-counting method of multifractal analysis of wavefunction on a Lattice
# Only works for cubic (L X L X L) lattice and cubic box (b X b X b)

module MFA_BC

import Lattice.Lattice1D, Lattice.Lattice2D, Lattice.Lattice3D,
    Lattice.index, Lattice.site

using LinearAlgebra
using Lattice
using PN

    @doc """
    Parameters of the box-counting method (only works for cubic lattice and cubic boxes).
    """
    struct BC_Params{T <: Integer, F <: AbstractFloat, L <: Lattice3D}
        lattice::L
        b::AbstractArray{T}
        q::AbstractArray{F}
        box_num::AbstractArray{T}
        box_size::AbstractArray{T}

        function BC_Params(ltc::L, b::AbstractArray{T}, q::AbstractArray{F}) where {T <: Integer, F <: AbstractFloat, L <: Lattice3D}
            if (ltc.L == ltc.M) && (ltc.M == ltc.N) && all(a -> ltc.L%a == 0, b)
                box_num = (ltc.L .÷ b).^3
                box_size = b.^3*ltc.U
                new{T, F, L}(ltc::L, vec(b)::AbstractArray{T}, vec(q)::AbstractArray{F}, box_num::AbstractArray{T}, box_size::AbstractArray{T})
            else
                error("Initialization failed. BC_Params only allow commensurate cubic boxes of cubic (L x L x L) lattice")
            end
        end
    end

    @doc """
    Indices of boxes of smaller cubic lattices of size (b X b X b)
    of cubic lattice of size (L X L X L), used for box-counting method.
    """
    function box_indices(ltc::Lattice3D, b::Int64)
        box_size = b^3 * ltc.U
        box_num = (ltc.L ÷ b)^3
        box_idx = Array{Int64}(undef, box_size, box_num)
        j = 1
        for l in 1:b:(ltc.L - b + 1), m in 1:b:(ltc.L - b + 1), n in 1:b:(ltc.L - b + 1)
            idx = Int64[] # indices of jth box
            for lj in 1:b, mj in 1:b, nj in 1:b, uj in 1:ltc.U
                box_site = ((l+lj-1), (m+mj-1), (n+nj-1), uj)
                push!(idx, index(ltc, box_site))
            end
            box_idx[:,j] = idx
            j += 1
        end
        return box_idx
    end

    @doc """
    Box count method for scaling analysis of GIPR. Compute scaling of GIPR of (L X L X L) cubic lattice using Box-counting method
    (TODO: Generalize to L X M X N boxes.)
    """
    function gipr(ψ, p::BC_Params)
        b = p.b; q = p.q
        out = Array{Float64}(undef,  length(b), length(q))
        for (i, bi) in enumerate(b)
            box_idx = box_indices(p.lattice, bi)
            ψ_boxed = Vector{Float64}(undef, p.box_num[i])
            for j in 1:p.box_num[i]
                ψj = @view ψ[box_idx[:,j]] # jth box
                ψ_boxed[j] =  sum(x -> abs2(x), ψj)
            end

            for k in 1:length(q)
                out[i, k] = sum(x -> x^q[k], ψ_boxed)
            end
        end
        return out
    end

    @doc """
    Fit the exponents of GIPR.
    """
    function fit_exp_gipr(gipr::AbstractArray{Float64, 2}, p::BC_Params; y_intercept = false)
        @assert size(gipr, 1) == length(p.b) && size(gipr, 2) == length(p.q)
        X = [ones(length(p.b)) log10.(p.lattice.L ./ p.b)]
        Y = log10.(gipr)

        v = X\Y
        if y_intercept
            return v
        else
            return v[2,:]
        end
    end

    @doc """
    legendre transform of (q,τ(q)).
    """
    function mf_spectrum(τ_q, p::BC_Params)
        τ_q = @view τ_q[end,:]
        α = diff(τ_q)./(p.q[2] - p.q[1])
        f_α = [p.q[i]*α[i] - τ_q[i] for i in 1:length(p.q)-1]

        return α, f_α
    end


    @doc """
    Parameters of the box-counting method (only works for cubic lattice and cubic boxes).
    """
    struct BC_Params2D{T <: Integer, F <: AbstractFloat, L <: Lattice2D}
        lattice::L
        b::AbstractArray{T}
        q::AbstractArray{F}
        box_num::AbstractArray{T}
        box_size::AbstractArray{T}

        function BC_Params2D(ltc::L, b::AbstractArray{T}, q::AbstractArray{F}) where {T <: Integer, F <: AbstractFloat, L <: Lattice2D}
            if (ltc.M == ltc.N) && all(a -> ltc.M%a == 0, b)
                box_num = (ltc.M .÷ b).^2
                box_size = b.^2*ltc.U
                new{T, F, L}(ltc::L, vec(b)::AbstractArray{T}, vec(q)::AbstractArray{F}, box_num::AbstractArray{T}, box_size::AbstractArray{T})
            else
                error("Initialization failed. BC_Params only allow commensurate cubic boxes of cubic (L x L x L) lattice")
            end
        end
    end

    @doc """
    Indices of boxes of smaller cubic lattices of size (b X b X b)
    of cubic lattice of size (L X L X L), used for box-counting method.
    """
    function box_indices(ltc::Lattice2D, b::Int64)
        box_size = b^2 * ltc.U
        box_num = (ltc.L ÷ b)^2
        box_idx = Array{Int64}(undef, box_size, box_num)
        j = 1
        for m in 1:b:(ltc.L - b + 1), n in 1:b:(ltc.L - b + 1)
            idx = Int64[] # indices of jth box
            for mj in 1:b, nj in 1:b, uj in 1:ltc.U
                box_site = ((m+mj-1), (n+nj-1), uj)
                push!(idx, index(ltc, box_site))
            end
            box_idx[:,j] = idx
            j += 1
        end
        return box_idx
    end
    export BC_Params, box_indices, gipr, fit_exp_gipr, mf_spectrum

end
