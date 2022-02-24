module Lattice
    using LinearAlgebra
    using SparseArrays

    @doc """
    1D Lattice struct
    """
    struct Lattice1D{T <: Int64, F <: Float64}
        N::T
        U::T
        index::Array{T, 2}
        sites::Array{NTuple{2, T}}
        H::SparseMatrixCSC{Complex{F}, Int64}

        function Lattice1D(N, U)
            index = Array{Int64}(undef, N, U)
            i = 1
            for n in 1:N, u in 1:U
                index[n,u] = i
                i += 1
            end
            sites = NTuple{2, Int64}[]
            for n in 1:N, u in 1:U
                push!(sites, (n, u))
            end

            H = spzeros(Complex{Float64}, N*U, N*U)
            new{Int64, Float64}(N, U, index, sites, H)
        end
    end

    @doc """
    2D Lattice struct
    """
    struct Lattice2D{T <: Int64, F <: Float64}
        M::T
        N::T
        U::T
        index::Array{T, 3}
        sites::Array{NTuple{3, T}}
        H::SparseMatrixCSC{Complex{F}, Int64}

        function Lattice2D(M, N, U)
            index = Array{Int64}(undef, M, N, U)
            i = 1
            for m in 1:M, n in 1:N, u in 1:U
                index[m,n,u] = i
                i += 1
            end
            sites = NTuple{3, Int64}[]
            for i in 1:M, j in 1:N, k in 1:U
                push!(sites, (i, j, k))
            end

            H = spzeros(Complex{Float64}, M*N*U, M*N*U)
            new{Int64, Float64}(M, N, U, index, sites, H)
        end
    end
    @doc """
    2D Lattice struct
    """
    struct Lattice3D{T <: Int64, F <: Float64}
        L::T
        M::T
        N::T
        U::T
        index::Array{T, 4}
        sites::Array{NTuple{4, T}}
        H::SparseMatrixCSC{Complex{F}, Int64}

        function Lattice3D(L, M, N, U)
            index = Array{Int64}(undef, L, M, N, U)
            i = 1
            for l in 1:L, m in 1:M, n in 1:N, u in 1:U
                index[l,m,n,u] = i
                i += 1
            end
            sites = NTuple{4, Int64}[]
            for l in 1:L, m in 1:M, n in 1:N, u in 1:U
                push!(sites, (l,m,n,u))
            end

            H = spzeros(Complex{Float64}, L*M*N*U, L*M*N*U)
            new{Int64, Float64}(L, M, N, U, index, sites, H)
        end
    end
    #
    # function visualize(ltc::Lattice2D)
    #     vts = ltc.sites
    #     p = scatter(vts, markersize = 5)
    #
    #     for i in 1:length(vts), j in 1:length(vts)
    #         if round(abs(ltc.H[i,j]), digits = 12) > 0 && [vts[i][1] vts[j][1]] != [1 ltc.M] && [vts[i][1] vts[j][1]] != [ltc.M 1]  && ([vts[i][2] vts[j][2]] != [1 ltc.N]) && ([vts[i][2] vts[j][2]] != [ltc.N 1])
    #             plot!(p, [vts[i]; vts[j]], color = "black", legend = false)
    #             println("$(vts[i]),$(vts[j]): $(ltc.H[i,j])")
    #         end
    #     end
    #     display(p)
    # end


    export Lattice2D, Lattice3D
end
