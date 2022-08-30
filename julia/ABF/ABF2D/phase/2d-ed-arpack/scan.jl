using LinearAlgebra
using Arpack
using DataFrames, CSV
using ArgParse, JSON
using GIPR
using Lattices, PN
using Parameters
using Statistics
using LatticeModels
using Random, Distributions
include("./phase.jl")

function midpoints(edges)
    mids = Array{eltype(edges)}(undef, length(edges)-1)
    for i in 1:length(edges)-1
        mids[i] = edges[i] + edges[i+1]
    end

    return mids/2
end

function readconfig(d::Dict)
    E_max = d["E_max"]
    E_min = d["E_max"]
    R = d["R"]
    nev = d["nev"]
    L = d["L"]
    l = d["l"]
    V = d["W"]
    θ = d["th"]
    cutoff = d["cutoff"]
    bins = d["bins"]
    seed = d["seed"] 
    return ScanParameters(E_max=E_max, E_min=E_min, R=R, nev=nev, L=L, l=l, V=V, θ=θ, cutoff=cutoff, bins=bins,seed=seed)
end


@with_kw struct ScanParameters
    E_max::Float64 = 0.3
    E_min::Float64 = 0.1
    R::Int = 5 
    nev::Int = 10
    L::Int = 80 
    l::Int = 8 
    V::Float64 = 1.
    θ::Float64 = 0.25
    cutoff::Float64 = 0.9
    bins::Int = 21
    seed::Int = 1234
end

function scan_2d(p::ScanParameters)
    @unpack E_max, E_min, R, nev, L, l, V, θ, cutoff, bins, seed = p
    rng = MersenneTwister(seed)
    # set the energy bins
    E_max  = E_max*cutoff
    E_min  = E_min
    E_edges= range(E_min, E_max, length=bins+1)
    E_mids = midpoints(E_edges) 
    # create box-counting indices
    ltc = Lattice2D(L, L, 1)
    box_inds = box_indices(ltc, l)
    iprs_full=Array{Float64}(undef, R, bins)

    for r in 1:R
        H = Hermitian(ham_sf(L=L, V=V, θ=θ, rng=rng))
        for i in 1:bins
            E, psi=eigs(H, nev=nev, sigma=E_mids[i])
            E = real.(E)

            limiter = E_edges[i] .< E .< E_edges[i+1]
            psi = psi[:, findall(limiter)]
            p = abs2.(psi)            
            p_coarse = box_coarse(p, box_inds)
            iprs_full[r, i] = mean(compute_iprs(p_coarse, density=true))
        end
           
    end
    tilde_tau = log.(mean(iprs_full, dims=1))./log(l/L)
    return E_mids, tilde_tau 
end

function main(ARGS)
    opts = ArgParseSettings(description="Scan and compute pn for all parameters of nu=2 ABF")
    @add_arg_table! opts begin
    "c"
        help = "configuration"
        arg_type = AbstractString
        required = false
    end
    # Parse the arguments
    popts   = parse_args(opts)
    if popts["c"] != nothing
        popts   = parse_args(opts)
        config  = JSON.parsefile(popts["c"])
        p = readconfig(config)
    else
        p = ScanParameters() 
    end

    fn = "tau-L$(p.L)-th$(p.θ).csv"
    @time E, tau = scan_2d(p)
    CSV.write(fn, DataFrame(E = E, tau = vec(tau)))
end

main(ARGS)
