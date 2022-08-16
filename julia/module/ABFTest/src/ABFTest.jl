module ABFTest
    import Lattice.Lattice1D, Lattice.Lattice2D, Lattice.Lattice3D,
        Lattice.index, Lattice.site
    using LinearAlgebra
    using SparseArrays

    include("ABF1D.jl")
    include("ABF2D.jl")
    include("ABF3D.jl")

    function project!(H_p::AbstractArray,H::AbstractArray)
        H_p .= H[2:2:end, 2:2:end]
    end

    function project(H::AbstractArray)
        H_p = H[2:2:end, 2:2:end]
        return H_p
    end

    export project!, project, Lattice
end
