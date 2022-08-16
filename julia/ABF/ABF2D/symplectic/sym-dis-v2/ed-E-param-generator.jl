using LinearAlgebra, SparseArrays, KrylovKit
using LinearMaps
using Random
using DataFrames, CSV, DelimitedFiles
using ArgParse, JSON
# Custom modules
using ABFSym
using Lattice
using PN
using Glob
include("./ed-params.jl") # read parameters from configuration file
include("./ed-main.jl")

function estimate_bw(p, θ, W, L, rng)
    ltc = Lattice2D(L, L, 4)
    H, U = ham_fe(ltc, -2, 0, θ) # Fully entangled hamiltonian
    H = convert.(ComplexF64, H)
    H_dis = makesym2d(ltc, H, p.V1, p.V2, rng = rng)
    D = Diagonal(dis(L^2, W, rng))
    @views H_prj = project(U'*(H_dis + D)*U)
    vals, psi, info = eigsolve(Hermitian(H_prj), size(H_prj, 1), 1, :LM, ishermitian = true, krylovdim = 30)
    return 2*abs(maximum(vals))
end

# @doc"""
# p: Params
# θ: angle
# W: disorder strength
# L: system size for estimating bandwidth.
# E_crop: crop the end of the spectrum. Recommend 0.9
# E_bin_width: portion given by 1/E_bin_width of the small energy window around energy center.
# rng: random number generator
# """
function auto_energy_params(p, θ, W, L, E_crop, E_bin_width, rng)
    BW = estimate_bw(p, θ, W, L, rng)
    E_c = range(0.0001, BW/2*E_crop, length = length(p.E))
    E_del = (E_c[2] - E_c[1])/E_bin_width
    return BW, E_c, E_del
end

function energy_param_generator(p, θ, W, L, E_crop, E_bin_width, rng)
    if p.bw_auto
        BW, E_c, E_del = auto_energy_params(p, θ, W, L, E_crop, E_bin_width, rng)
    else
        if length(p.E) != 1
            BW = p.E[end] - p.E[1]
            E_c = p.E
            E_del = (E_c[2] - E_c[1])/p.E_bin_width
        elseif length(p.E) == 1
            BW = estimate_bw(p, θ, W, L, rng)
            E_c = p.E[1] + 1E-12
            E_del = BW/p.E_bin_width
        end
    end
    return BW, E_c, E_del
end

function abf2d_energy_param_scan(p::Params)
    rng = MersenneTwister(p.seed)
    BW = Array{Float64}(undef, length(p.θ), length(p.W))
    E_c = Array{Float64}(undef, length(p.E), length(p.θ) , length(p.W))
    E_del =   
    for j in 1:length(p.θ), jj in 1:length(p.W)
        BW, E_c, E_del = energy_param_generator(p, p.θ[j], p.W[jj], 200, 0.90, 4, rng)
    end
end
