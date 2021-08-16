using Lattice
using SparseArrays

@doc """
Nearest neighbor strip with periodic boundary conditions imposed on both the width and length.

ltc: 2D lattice struct Lattice2D(M, N, U). by default Lattice2D(3, 10, 1).
(M: width of the strip. M > 2, N:length of the strip. N > 2, U = 1).

t: hopping strength default is 1..
"""
function ham_strip(;M = 3, N = 10, t = 1)
    @assert M > 2 || N > 2
    ltc = Lattice2D(M, N, 1)

    I = Int64[]; J = Int64[];
    for m in 1:ltc.M, n in 1:ltc.N
        push!(I, index(ltc, (m, n, 1)))
        push!(J, index(ltc, (m, n+1, 1)))
        push!(I, index(ltc, (m, n+1, 1)))
        push!(J, index(ltc, (m, n, 1)))

        push!(I, index(ltc, (m, n, 1)))
        push!(J, index(ltc, (m+1, n, 1)))
        push!(I, index(ltc, (m+1, n, 1)))
        push!(J, index(ltc, (m, n, 1)))
    end

    V = t*ones(Float64, length(I))
    H = sparse(I, J, V, ltc.N*ltc.M, ltc.N*ltc.M)
    return H
end

@doc """
Nearest neighbor strip with periodic boundary conditions imposed on both the width and length.

ltc: 2D lattice struct Lattice2D(M, N, U). by default Lattice2D(3, 10, 1).
(M: width of the strip. M > 2, N:length of the strip. N > 2, U = 1).

t: hopping strength default is 1..
"""
function ham_chain(;N = 10, t = 1)
    ltc = Lattice1D(N, 1)

    I = Int64[]; J = Int64[];
    for n in 1:ltc.N
        push!(I, index(ltc, (n, 1)))
        push!(J, index(ltc, (n+1, 1)))
        push!(I, index(ltc, (n+1, 1)))
        push!(J, index(ltc, (n, 1)))
    end

    V = t*ones(Float64, length(I))
    H = sparse(I, J, V, ltc.N, ltc.N)
    return H
end


@doc """
Nearest neighbor bar with periodic boundary conditions on width and length.

ltc: 3D lattice struct Lattice2D(L, M, N, U). by default Lattice3D(3, 3, 10, 1).
(L, M: width of the bars. L == M, N:length of the strip. N > 2, U = 1).

t: hopping strength default is 1..
"""
function ham_bar(;M = 3, N = 10, t = 1)
    ltc = Lattice3D(M, M, N, 1)
    @assert ltc.L > 2 || ltc.M > 2 || ltc.N > 2
    @assert ltc.L == ltc.M
    I = Int64[]; J = Int64[];
    for l in 1:ltc.L, m in 1:ltc.M, n in 1:ltc.N
        push!(I, index(ltc, (l, m, n, 1)))
        push!(J, index(ltc, (l, m, n+1, 1)))
        push!(I, index(ltc, (l, m, n+1, 1)))
        push!(J, index(ltc, (l, m, n, 1)))

        push!(I, index(ltc, (l, m, n, 1)))
        push!(J, index(ltc, (l, m+1, n, 1)))
        push!(I, index(ltc, (l, m+1, n, 1)))
        push!(J, index(ltc, (l, m, n, 1)))

        push!(I, index(ltc, (l, m, n, 1)))
        push!(J, index(ltc, (l+1, m, n, 1)))
        push!(I, index(ltc, (l+1, m, n, 1)))
        push!(J, index(ltc, (l, m, n, 1)))
    end

    V = t*ones(Float64, length(I))
    H = sparse(I, J, V, ltc.L*ltc.N*ltc.M, ltc.L*ltc.N*ltc.M)
    return H
end
