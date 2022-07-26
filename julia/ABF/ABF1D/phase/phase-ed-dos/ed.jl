using LinearAlgebra, SparseArrays, KrylovKit
using Statistics, StatsBase
using LinearMaps
using Random
using DataFrames, CSV
# Custom modules
using ABF
using Lattice
using PN
using JLD
using Binning
include("./params.jl") # read parameters from configuration file
include("./boxcount.jl")

# Functions related to lanczos matrix diagonalization
function construct_linear_map(A)
    F = factorize(A)
    LinearMap{eltype(A)}((y, x) -> ldiv2!(y, F, x), size(A, 1), ismutating = true, ishermitian = true)
end
function ldiv2!(y, F, x)
    y .= F\x
end

function phase_dis(H; V = 1., rng = nothing)
    D = convert.(ComplexF64, H)
    if rng == nothing
        rng = MersenneTwister()
    end
    rows = rowvals(D)
    vals = nonzeros(D)
    m, n = size(D)
    for j = 1:n
       for i in nzrange(D, j)
          row = rows[i]
          # println("$row ,", "$j")
          if row > j
              vals[i] = im*V*vals[i]*(rand(rng) .- 0.5)
          elseif row <= j
              vals[i] = 0.
          end
       end
    end
    return D + D'
end

function auto_estimate_bw(θ::Float64, V::Float64, W::Float64, rng::AbstractRNG)
    L = 5000

    H, U = ham_fe(Lattice1D(L, 2), -2., 0., θ) # Fully entangled hamiltonian
    H = convert.(ComplexF64, H)
    max_E = Float64[]

    for i in 1:20
        D_phase = phase_dis(H, V = 1, rng = rng)
        D_onsite = Diagonal(rand(rng, 2L) .- 0.5)
        @views H_prj = project(U'*(V*D_phase + W*D_onsite)*U)
        vals, _, _= eigsolve(Hermitian(H_prj), size(H_prj, 1), 1, :LM, ishermitian = true, krylovdim = 30)
        push!(max_E, maximum(vals))
    end    
    return mean(max_E)
end

function generate_energy_file(p::Params)
    mkpath(p.energy_path)
    @unpack L, θ, W, seed, R, V, num_blas = p
    rng = MersenneTwister(seed)
    max_E = Array{Float64}(undef, length(θ), length(W)) 
    for i in 1:length(θ), j in 1:length(W)
        max_E[i, j] = auto_estimate_bw(θ[i], V, W[j], rng)
    end
    JLD.save(p.energy_path*"/max_E.jld", "θ",θ,"W", W,"max_E", max_E)
end

function get_energy_file(p::Params)
    done = false
    local dict
    while !done
        if isfile(p.energy_path*"/max_E.jld")
            dict = JLD.load(p.energy_path*"/max_E.jld")
            criteria = dict["θ"] == p.θ && dict["W"] == p.W
            if !criteria
                error("The parameters of the energy file and the current parameters does not match")
            end
            done = true 
        else
            println("Energy file is not found. Genereting energy files")
            generate_energy_file(p)
        end
    end
    return dict["max_E"]
end

function midpoints(e_edges)
    return [(e_edges[i] + e_edges[i+1])/2 for i in 1:length(e_edges)-1]
end

function abf_scan(p::Params)
    @unpack L, θ, W, E_edges, seed,R, V, num_blas, energy_path = p
    mkpath("L$(L)")
    rng = MersenneTwister(seed) 
    ltc = Lattice1D(L, 2)
    max_E = get_energy_file(p)
    println("Energ file is obtained")
    for i in 1:length(θ)
        H,U=ham_fe(ltc, -2., 0., θ[i]) # Fully entangled hamiltonian
        H=convert.(ComplexF64, H)
        for j in 1:length(W)
            E_edges=max_E[i,j]*p.E_edges
            E_midpoints=midpoints(E_edges)
            df = DataFrame(E=E_midpoints,E_idx=collect(1:length(E_midpoints)),weight=zeros(Int64,length(E_midpoints)),dos=zeros(Float64,length(E_midpoints)))
             
            for r in 1:R # Realizations
                D_onsite = Diagonal(rand(rng, 2p.L) .- 0.5)
                D_phase = phase_dis(H, V = 1, rng = rng)
                @views H_prj = project(U'*(W[j]*D_onsite + V*D_phase)*U)
                droptol!(H_prj, 1E-14)
                e = eigvals!(Hermitian(Matrix(H_prj)))
                lbl = binning(e, E_edges)
                binned_histogram!(df.weight, e, E_edges)
            end
            df.dos = df.weight ./ diff(E_edges)
            savefn = "L$(L)/Th$(lpad(i, 2, "0"))_W$(lpad(j, 2, "0")).csv"
            CSV.write(savefn, df)
            println("Written:", savefn)
        end
    end
end

