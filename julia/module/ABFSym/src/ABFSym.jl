module ABFSym
    using Lattices
    using LinearAlgebra
    using SparseArrays

    include("ABF1D_Sym.jl")
    include("ABF2D_Sym.jl")
    include("ABF3D.jl")

    function project(H::AbstractArray)
        N = size(H, 1)÷4
        ind = Int64[]
        for i in 1:N
            push!(ind, 4(i-1) + 3)
            push!(ind, 4(i-1) + 4)
        end
        H_p = H[ind, ind]
        return H_p
    end

    export project
end
