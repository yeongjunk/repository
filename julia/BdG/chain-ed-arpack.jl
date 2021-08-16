using Lattice
using SparseArrays
using LinearAlgebra
using Random
using PN
using Arpack
using Statistics
using Parameters

function bdg_chain(;N = 3, J = 1., ga = 1., W = 1., rng = nothing)
    ltc = Lattice1D(N, 1)
    row = Int64[]; col = Int64[]; V = Float64[];
    if rng == nothing
        rng = MersenneTwister()
    end

    ε = W*(rand(rng, N) .- 0.5)
    for n in 1:N
        ζ_n = (1 - ε[mod1(n-1, N)]/(2ga) + 1 - ε[mod1(n+1, N)]/(2ga))/(1 - ε[n]/(2*ga))
#         println(ζ_n)
        c_S_n = J*2(ga-ε[n])*ζ_n
        c_S_n1 = -J*2(ga-ε[n])

        push!(row, index(ltc, (n, 1)))
        push!(col, index(ltc, (n, 1)))
        push!(V, c_S_n)

        push!(row, index(ltc, (n, 1)))
        push!(col, index(ltc, (n+1, 1)))
        push!(V, c_S_n1)

        push!(row, index(ltc, (n, 1)))
        push!(col, index(ltc, (n-1, 1)))
        push!(V, c_S_n1)
    end
    return sparse(row, col, V, N, N)
end

## scan parameter
R = 10; J = 1; N = 100000; nev = 1; ga = 2.01; W = 4
λ2 = (0.:0.1:6.).^2
##
@with_kw struct ScanParameters
    R::Int = 1
    J::Float64 = 1.
    N::Int = 10000
    nev::Int = 1
    ga::Float64 = 50
    W::Float64 = 4
    λ2::Array{Float64}
end

function scan_bdg(p::ScanParameters)
    @unpack R, J, N, nev, ga, W, λ2 = p
    pn = similar(λ2)
    E = similar(λ2)
    fill!(pn, 0)
    fill!(E, 0)
    for r in 1:R
        H = bdg_chain(N = N, W = W, ga = ga)
        for i in 1:length(λ2)
            λ_temp, psi = eigs(H, nev = nev, sigma = λ2[i])
            E[i] += λ_temp[1] # only takes nearest one
            pn[i] += mean(compute_pns(psi))
        end
    end
    E ./= R
    λ = sqrt.(abs.(E .+ 1E-11))
    pn ./= R
    return λ, pn
end
