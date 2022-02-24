# Useful functions for scanning

using Random
using LinearAlgebra
using StatsBase
include("abf2_ham.jl")

@doc """
A struct where parameters are stored. This is the input parameters for scanning.
θI and θII are given in unit of π
"""
struct Parameters{F<:AbstractFloat, N<:Int64}
    N::N # unit cell
    W_min::F
    W_max::F
    W_num::N
    θ::F
    ϕ_min::F
    ϕ_max::F
    ϕ_num::N
    seed::N
    R::N
    PN_bin_num::N
    PN_bin_min::N
    PN_bin_max::N
    E_bin_num::N
end

@doc """
Read configuration JSON file, and construct a Parameters struct from it.
"""
function readconfig(config::Dict)
    p = Parameters(
        config["N"],
        config["W_min"],
        config["W_max"],
        config["W_num"],
        config["theta"],
        config["phi_min"],
        config["phi_max"],
        config["phi_num"],
        config["seed"],
        config["R"],
        config["PN_bin_num"],
        config["PN_bin_min"],
        config["PN_bin_max"],
        config["E_bin_num"]
        )

    return p
end

function expand_params(p::Parameters)
    ϕ = collect(range(p.ϕ_min, p.ϕ_max, length = p.ϕ_num))
    W = collect(range(p.W_min, p.W_max, length = p.W_num))
    return W, ϕ
end

@doc """
computation of PN
"""
function compute_pn(eigvect::AbstractVector{ComplexF64})
    partn = sum(x->abs(x)^4, eigvect)
    return 1 / partn
end

function compute_pn_col!(dest::AbstractVector{Float64}, eigvect::AbstractArray{ComplexF64, 2})
    @assert length(dest) == size(eigvect, 2)
    for i in 1:length(dest)
        dest[i] = compute_pn(@view eigvect[:,i])
    end
end

function sum_bin(E::AbstractArray{Float64,1}, PN::AbstractArray{Float64,1}, min::Float64, max::Float64, bins::Int64)
    sum_bin = zeros(Float64, bins)
    N_bin = zeros(Int64, bins)
    step = (max-min)/bins
    for i = 1:length(E)
        j = convert(Int64, div(E[i]-min,step)) + 1
        if j <= bins
            sum_bin[j] += PN[i]
            N_bin[j] +=1
        end
    end
    return N_bin, sum_bin
end

function normalize_E(E::Float64,E0::Float64,W::Float64)
    E = (E - E0) / W
end

function revert_E(norm_E,W,E0)
     E = norm_E*W + E0
 end



@doc """
Scan histogram, mean, variance of PN of R, scanning W, θ.
"""
function scan_pn(p::Parameters)

    println("Initializing...")
    # setup parameters
    W, ϕ = expand_params(p)
    # Initialize output arrays
    println("Initializing output space...")
    num_samp    = zeros(Int64, p.ϕ_num, p.W_num, p.E_bin_num) # of samples for each W, θ
    ΣPN  = zeros(Float64, p.ϕ_num, p.W_num, p.E_bin_num)
    ΣPN² = zeros(Float64, p.ϕ_num, p.W_num, p.E_bin_num)

    println("Initializing histograms...")
    # initialize histogram array
    PN_bins = range(p.PN_bin_min, p.PN_bin_max, length = p.PN_bin_num) # Histogram bins
    E_bins = range(-0.5, 0.5, length = p.E_bin_num)
    PN_hist = Array{Histogram}(undef, p.ϕ_num, p.W_num)
    for i in eachindex(PN_hist)
        PN_hist[i] =  fit(Histogram, ([],[]), (E_bins, PN_bins))
    end

    H_det = Array{Float64}([1 0 0; 0 0 0; 0 0 -1])
    U = unitary_block(p.θ)
    T = uc_redef(p.N)
    U_full = blockdiagonal(U, p.N)
    x = 0
    println("Initialized.")

    println("Start scanning...")
    for n = 1:p.R #Iteration for disorder realizations
        rng = MersenneTwister(p.seed + n)
        r_arr = rand(rng, 3p.N) .- 0.5
        progress = n/p.R
        if progress > 0.1x
            println("$(n) realization out of $(p.R)")
            x += 1
        end

        Threads.@threads for i = 1:length(ϕ)

            # H_det1 = U*H_det*U'
            # H_det1[2,3] = H_det1[2,3]*exp(π*im*ϕ[i])
            # H_det1[3,2] = conj(H_det[2,3])
            t = exp(2π*im*ϕ[i])
            H_det_block = Array{ComplexF64}([0 -t -conj(t); -conj(t) 0 -t; -t -conj(t) 0])
            H = blockdiagonal(H_det_block, p.N)
            H = T*H*T'
            H = U_full*H*U_full'

            E0 = -2. * cos(2. * pi * ϕ[i])
            for j = 1:length(W) # Iteration for W
                r_onsite = W[j]*r_arr

                H_dis = H + Diagonal(r_onsite)
                eig = eigen(Hermitian(Array(H_dis)),E0-W[j]/2,E0+W[j]/2)
                eig_num = length(eig.values)
                PN = Vector{Float64}(undef, eig_num)
                compute_pn_col!(PN, eig.vectors)

                E = eig.values
                E = normalize_E.(E, E0, W[j]) # normalized energy

                PN_hist_temp = fit(Histogram, (E, PN), (E_bins, PN_bins))
                num_temp, ΣPN_temp = sum_bin(E, PN, -0.5, 0.5, p.E_bin_num)
                _, ΣPN²_temp = sum_bin(E, PN.^2, -0.5, 0.5, p.E_bin_num)

                num_samp[i, j, :] .+= num_temp
                ΣPN[i, j,:] .+= ΣPN_temp
                ΣPN²[i, j, :] .+= ΣPN²_temp

                merge!(PN_hist[i,j], PN_hist_temp)
            end
        end
    end
    PN_mean = ΣPN ./ num_samp
    PN²_mean = ΣPN² ./ num_samp
    PN_var = num_samp ./ (num_samp .-1) .* (PN²_mean - PN_mean.^2) #

    return num_samp, PN_mean, PN_var, PN_hist
end
