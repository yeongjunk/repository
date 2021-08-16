# Useful functions for scanning

using Random
using LinearAlgebra
using StatsBase
include("abf2_ham.jl")
include("abf2_block.jl")

@doc """
A struct where parameters are stored. This is the input parameters for scanning.
θI and θII are given in unit of π
"""
struct Parameters{F<:AbstractFloat,N<:Int64}
    E::F #energy
    E_range::F #Should be given as ratio (E_range = 0.1 means +-10% of entire region)
    N::N # unit cell

    W_min::F
    W_max::F
    W_num::N

    θI_min::F
    θI_max::F
    θI_num::N

    θII_min::F
    θII_max::F
    θII_num::N

    seed::N
    R::N #realizations
    bin_num::N
    sample_max::N # prevention from calculating too many samples at weak disorder
end

@doc """
Read configuration JSON file, and construct a Parameters struct from it.
"""
function readconfig(config::Dict)
    p = Parameters(
        config["E"],
        config["E_range"],
        config["N"],
        config["W_min"],
        config["W_max"],
        config["W_num"],
        config["theta_min"][1],
        config["theta_max"][1],
        config["theta_num"][1],
        config["theta_min"][2],
        config["theta_max"][2],
        config["theta_num"][2],
        config["seed"],
        config["R"],
        config["sample_max"],
        config["bin_num"],
        config["bin_min"],
        config["bin_max"]
    )

    return p
end

function expand_params(p::Parameters)
    θI = collect(range(p.θI_min, p.θI_max, length = p.θI_num))
    θII = collect(range(p.θII_min, p.θII_max, length = p.θII_num))
    W = collect(range(p.W_min, p.W_max, length = p.W_num))

    return W, θI, θII
end

@doc """
computation of PN
"""
function compute_pn(eigvect::AbstractVector{Float64})
    partn = sum(x->abs(x)^4, eigvect)
    return 1 / partn
end

function compute_pn_col!(dest::AbstractVector{Float64}, eigvect::AbstractArray{Float64, 2})
    @assert length(dest) == size(eigvect, 2)
    for i in 1:length(dest)
        dest[i] = compute_pn(@view eigvect[:,i])
    end
end

@doc """
Scan histogram, mean, variance of PN of R, scanning W, θ.
"""
function scan_pn(p::Parameters; logscale = false)
    # setup parameters
    W, θ, _ = expand_params(p)

    if logscale
        W = 10 .^W
        println("logscale is set.")
    end

    # Initialize output arrays
    smp_num = zeros(Int64, p.θI_num, p.W_num) # of samples for each W, θ
    pn_sum= zeros(Float64, p.θI_num, p.W_num)
    pnsq_sum = zeros(Float64, p.θI_num, p.W_num)
    # initialize histogram array
    pn_hist = Array{Histogram}(undef, p.θI_num, p.W_num)
    bins = range(p.bin_min, p.bin_max, length = p.bin_num) # Histogram bins
    for i in eachindex(pn_hist)
        pn_hist[i] =  fit(Histogram,[], bins)
    end

    # allocate space for hamiltonain matrices
    td_num = Threads.nthreads()
    H = zeros(Float64, 2p.N, 2p.N, td_num) # individual hamiltonian matrix storage for each threads.

    for n = 1:p.R #Iteration for each disorder realizations
        rng = MersenneTwister(p.seed + n)
        r_arr = rand(rng, 2p.N) .- 0.5
        for i = 1:length(θ)
            for t = 1:td_num
                overwrite_ham_abf2!(view(H, :, :, t), θ[i], θ[i], p.N)
            end
            Threads.@threads for j = 1:length(W) # Iteration for W
                if smp_num[i,j] > p.sample_max
                    continue
                end

                r_onsite = W[j]*r_arr
                H_dis = view(H, :, :, Threads.threadid())
                add_diag!(H_dis, r_onsite)
                eig = eigen(Symmetric(H_dis), p.E - 2p.E_range, p.E + 2p.E_range)
                sub_diag!(H_dis, r_onsite) #remove disorder for recycling

                eig_num = length(eig.values)
                if eig_num == 0 # there is no sample within the energy window
                    continue
                end

                pn = Vector{Float64}(undef, eig_num)
                compute_pn_col!(pn, eig.vectors)

                temp_pn_hist = fit(Histogram, pn, bins)

                smp_num[i,j] += eig_num
                pn_sum[i, j] += sum(pn)
                pnsq_sum[i,j] += sum(x -> x^2, pn)
                merge!(pn_hist[i,j], temp_pn_hist)

            end
        end
    end
    pn_mean = pn_sum ./ smp_num
    pnsq_mean = pnsq_sum ./ smp_num
    pn_var = smp_num ./ (smp_num .-1) .* (pnsq_mean - pn_mean.^2) #
    println("mincount = $(findmin(smp_num)), maxcount = $(findmax(smp_num))")

    return smp_num, pn_mean, pn_var, pn_hist
end
