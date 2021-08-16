SAVEDIR = "/Users/pcs/data/cubic/"
using CSV, DataFrames
using Random
using LinearAlgebra
using SparseArrays

# Loadc custom modules
using Lattice
using PN

function cubic(L, t)
    M = L; N = L
    ltc = Lattice3D(L, M, N, 1)
    for l in 1:L, m in 1:M, n in 1:N
        ltc.H[ltc.index[mod1(l+1, L),m,n, 1],ltc.index[l,m,n, 1]] = -t
        ltc.H[ltc.index[mod1(l-1, L),m,n, 1],ltc.index[l,m,n, 1]] = -t
        ltc.H[ltc.index[l,mod1(m+1, M),n, 1],ltc.index[l,m,n, 1]] = -t
        ltc.H[ltc.index[l,mod1(m-1, M),n, 1],ltc.index[l,m,n, 1]] = -t
        ltc.H[ltc.index[l,m,mod1(n+1, N), 1],ltc.index[l,m,n, 1]] = -t
        ltc.H[ltc.index[l,m,mod1(n-1, N), 1],ltc.index[l,m,n, 1]] = -t
    end
    return ltc
end

q = range(-3, 3, length = 31)
L = 5:5:25
R = [1000 ;400; 100; 50; 15]

W = 10.
t = 1.

rng = MersenneTwister(12513)


q_str = string.(q)
L_str = string.(L)
W_str = replace(string(W), '.'=>"-")
t_str = replace(string(t), '.'=>"-")

for i in 1:length(L)
    l = cubic(L[i], t)
    df = DataFrame()
    for r in 1:R[i]
        D = spdiagm(0 => W*(rand(rng, L[i]^3) .- 0.5))
        l.H .+= D
        E, ψ = eigen(Hermitian(Matrix(l.H)))
        l.H .-= D

        gipr = [compute_iprs(ψ, q = q[j]) for j in 1:length(q)]

        df_t = DataFrame()
        df_t.E = E

        for j in 1:length(q)
            insertcols!(df_t, q_str[j] => gipr[j])
        end
        df_t.R = fill(r, length(E))

        append!(df, df_t)
    end
    fn = "L"*L_str[i]*"_W"*W_str*"_t"*t_str*".csv"
    CSV.write(SAVEDIR*fn, df)
    println("SAVED $(fn)")
end
