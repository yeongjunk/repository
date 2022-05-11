# A convenient tool for indexing vectors and matrices of lattices. 

module Lattice
    @doc """
    1D Lattice
    """
    struct Lattice1D{T<:Integer}
        N::T
        U::T
    end

    @doc """
    2D Lattice
    """
    struct Lattice2D{T<:Integer}
        M::T
        N::T
        U::T
    end

    @doc """
    3D Lattice
    """
    struct Lattice3D{T<:Integer}
        L::T
        M::T
        N::T
        U::T
    end


    @doc """
    Convert a lattice index to the lattice site (n,u)
    """
    function index(ltc::Lattice1D, (n,u)::NTuple{2, Integer})
        @assert u <= ltc.U
        N = ltc.N
        U = ltc.U

        n = mod1(n, N)
        i = (n - 1)*U + (u - 1) + 1
        return i
    end

    @doc """
    Convert lattice site (n,u) to lattice index
    """
    function site(ltc::Lattice1D, i::Integer)
        N = ltc.N
        U = ltc.U
        j = i-1
        # Compute zero-based index site
        n = div(j, U)
        u = mod(j, U)

        return (n + 1, u + 1) # convert to one-based site
    end


    @doc """
    Corresponding index of lattice site (m, n, u)
    """
    function index(ltc::Lattice2D, (m,n,u)::NTuple{3, Integer})
        @assert u <= ltc.U
        M = ltc.M
        N = ltc.N
        U = ltc.U

        m = mod1(m, M)
        n = mod1(n, N)
        i = (m - 1)*N*U + (n - 1)*U + (u - 1) + 1
        return i
    end

    @doc """
    lattice position (nmu) to matrix index
    """
    function site(ltc::Lattice2D, i::Integer)
        M = ltc.M
        N = ltc.N
        U = ltc.U
        j = i -1
        # Compute zero-based site
        m = div(j, N*U)
        n = div(j - m*N*U, U)
        u = mod(j, U)
        return (m + 1, n + 1, u + 1) #One-based site
    end

    @doc """
    Corresponding index of lattice site (m, n, u)
    """
    @inline function index(ltc::Lattice3D, (l, m, n, u)::NTuple{4, Integer})
        @assert u <= ltc.U
        L = ltc.L
        M = ltc.M
        N = ltc.N
        U = ltc.U

        l = mod1(l, L)
        m = mod1(m, M)
        n = mod1(n, N)
        i = (l - 1)*M*N*U + (m - 1)*N*U + (n - 1)*U + (u - 1) + 1
        return i
    end

    @doc """
    lattice position (l, m, n, u) to matrix index
    """
    function site(ltc::Lattice3D, i::Integer)
        j = i-1 # zero-base index
        L = ltc.L
        M = ltc.M
        N = ltc.N
        U = ltc.U
        # zero-base site
        l = div(j, M*N*U)
        m = div(j - l*M*N*U, N*U)
        n = div(j - l*M*N*U - m*N*U, U)
        u = mod(j,U)

        return (l+1, m+1, n+1, u+1)
    end

    export Lattice3D, Lattice2D, Lattice1D, index, site
end
