SAVEDIR = "/Users/pcs/data/square/"
using CSV, DataFrames
using Random
using LinearAlgebra
using SparseArrays

# Loadc custom modules
using Lattice
using PN

function square(L, t)
    ltc = Lattice2D(L, L, 1)

    num_sites = ltc.M*ltc.N*ltc.U
    I = Int64[]; J = Int64[]; V = Float64[]
    for m in 1:ltc.M, n in 1:ltc.N
        push!(I, index(ltc, (m, n, 1)))
        push!(J, index(ltc, (m, mod1(n-1, ltc.N),1)))
        push!(V, -t)

        push!(I, index(ltc, (m, n, 1)))
        push!(J, index(ltc, (m, mod1(n+1, ltc.N),1)))
        push!(V, -t)

        push!(I, index(ltc, (m, n, 1)))
        push!(J, index(ltc, (mod1(m-1, ltc.M),n,1)))
        push!(V, -t)

        push!(I, index(ltc, (m, n, 1)))
        push!(J, index(ltc, (mod1(m+1, ltc.M),n,1)))
        push!(V, -t)
    end
    return sparse(I, J, V, num_sites, num_sites)
end

q = range(0, 3, length = 31)
# q = range(2, 2, length = 1)

L = [20 60 100 140]
R = [3000 300 80 50]
# L = [20]
# R = [2]


W = range(1.5, 4.5, step = 0.5)
t = 1.

rng = MersenneTwister(12561)


q_str = string.(q)
L_str = string.(L)

@time for i in 1:length(L)
    H = square(L[i], t)
    df = DataFrame()
    for r in 1:R[i]
        for Wi in W
            D = spdiagm(0 => Wi*(rand(rng, L[i]^2) .- 0.5))
            H .+= D
            E, ψ = eigen(Hermitian(Matrix(H)), -0.1, 0.1)
            H .-= D

            gipr = [round.(compute_iprs(ψ, q = q[j]), sigdigits = 7) for j in 1:length(q)]
            df_t = DataFrame()
            df_t.E = round.(E, sigdigits = 9)

            for j in 1:length(q)
                insertcols!(df_t, q_str[j] => gipr[j])
            end
            df_t.R = fill(r, length(E))
            df_t.W = fill(Wi, length(E))
            append!(df, df_t)
        end
    end
    fn = "L"*L_str[i]*".csv"
    CSV.write(SAVEDIR*fn, df)
    println("SAVED $(fn)")
end
