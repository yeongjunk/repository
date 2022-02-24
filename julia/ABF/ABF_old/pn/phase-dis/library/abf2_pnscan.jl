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
struct Parameters{F<:AbstractFloat,N<:Int64, B<:Bool}
    E::F #energy
    E_window::F
    N::N # unit cell

    V_min::F
    V_max::F
    V_num::N
    V_log::B

    W::F

    θI_min::F
    θI_max::F
    θI_num::N

    θII_min::F
    θII_max::F
    θII_num::N

    # Statistics parameters
    seed::N
    R::N
    PN_bin_num::N
    PN_min::N
    PN_max::N
    E_bin_num::N
end

@doc """
Read configuration JSON file, and construct a Parameters struct from it.
"""
function readconfig(config::Dict)
    p = Parameters(
        config["E"],
        config["E_window"],
        config["N"],
        config["V_min"],
        config["V_max"],
        config["V_num"],
        config["V_log"],
        config["W"],
        config["theta_min"][1],
        config["theta_max"][1],
        config["theta_num"][1],
        config["theta_min"][2],
        config["theta_max"][2],
        config["theta_num"][2],
        config["seed"],
        config["R"],
        config["PN_bin_num"],
        config["PN_min"],
        config["PN_max"],
        config["E_bin_num"],
        )

    return p
end

function expand_params(p::Parameters)
    θI = collect(range(p.θI_min, p.θI_max, length = p.θI_num))
    θII = collect(range(p.θII_min, p.θII_max, length = p.θII_num))
    V = collect(range(p.V_min, p.V_max, length = p.V_num))

    return V, θI, θII
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

function compute_pn(H::AbstractArray{ComplexF64, 2}, E_min::Float64, E_max::Float64)
    eig = eigen(Hermitian(H), E_min, E_max)
    eig_num = length(eig.values)
    PN = Vector{Float64}(undef, eig_num)
    compute_pn_col!(PN, eig.vectors)
    return eig.values, PN
end


function sum_bin(E::AbstractVector{Float64}, PN::AbstractVector{Float64}, min::Float64, max::Float64, bins::Int64)
    sum_bin = zeros(Float64, bins)
    N_bin = zeros(Int64, bins)
    step = (max-min) / bins
    for i = 1:length(E)
        j = convert(Int64, div(E[i] - min,step) + 1)
        if 0 < j <= bins
            sum_bin[j] += PN[i]
            N_bin[j] += 1
        end
    end
    return N_bin, sum_bin
end

function sum_bin2(E::AbstractVector{Float64}, PN::AbstractVector{Float64}, min::Float64, max::Float64, bins::Int64)
    sum_bin = zeros(Float64, bins)
    step = (max-min) / bins
    for i = 1:length(E)
        j = convert(Int64, div(E[i] - min,step) + 1)
        if 0 < j <= bins
            sum_bin[j] += PN[i]
        end
    end
    return sum_bin
end


@doc """
Scan histogram, mean, variance of PN of R, scanning W, θ.
"""
function scan_pn_diag(p::Parameters; logscale::Bool)
    # setup parameters
    V, θI, _ = expand_params(p)
    if logscale
        V = 10 .^V
        println("logscale set...")
    end

    # Initialize output arrays
    num_samp    = zeros(Int64, p.θI_num, p.V_num, p.E_bin_num) # of samples for each W, θ
    ΣPN  = zeros(Float64, p.θI_num, p.V_num, p.E_bin_num)
    ΣPN² = zeros(Float64, p.θI_num, p.V_num, p.E_bin_num)

    # initialize histogram array
    PN_hist = Array{Histogram}(undef, p.θI_num, p.V_num)
    E_bins = range(-0.5, 0.5, length = p.E_bin_num)

    println("Output array initialized.")
    PN_max_scale = (p.PN_max/p.θI_max)*θI .+ 0.1
    for i in 1:p.θI_num
        for j in 1:p.V_num
            PN_bins = range(p.PN_min, PN_max_scale[i], length = p.PN_bin_num)
            PN_hist[i,j] =  fit(Histogram, ([],[]), (E_bins, PN_bins))
        end
    end

    H = zeros(Float64, 2p.N, 2p.N)
    for n = 1:p.R #Iteration for each disorder realizations
        rng = MersenneTwister(p.seed + n)
        r_on = p.W*(rand(rng, 2p.N) .- 0.5)
        r_hop_norm = rand(rng, 2p.N, 3) .- 0.5

        @Threads.threads for i = 1:length(θI)
            local H = create_ham_abf2(θI[i], θI[i], p.N)
            for j = 1:length(V) # Iteration for W
                r_hop = V[j]*r_hop_norm
                H_dis = add_diag(H, r_on)
                H_dis = convert(Array{ComplexF64,2}, H_dis)
                add_hop_dis!(H_dis, r_hop)

                (E, PN) = compute_pn(H_dis, -2., 0.)
                println("V: $(V[j]), θ: $(θI[i])"*", PN:$(minimum(PN))")
                E = (E .+ 1) / p.W
                PN_bins = range(p.PN_min, PN_max_scale[i], length = p.PN_bin_num)
                PN_hist_temp = fit(Histogram, (E, PN), (E_bins, PN_bins))
                num_samp_temp, ΣPN_temp = sum_bin(E, PN, -0.5, 0.5, p.E_bin_num)
                ΣPN²_temp = sum_bin2(E, PN.^2, -0.5, 0.5, p.E_bin_num)

                num_samp[i, j, :] .+= num_samp_temp
                ΣPN[i, j,:] .+= ΣPN_temp
                ΣPN²[i, j, :] .+= ΣPN²_temp

                merge!(PN_hist[i,j], PN_hist_temp)
            end
        end
    end

    PN_mean = ΣPN ./ num_samp
    PN²_mean = ΣPN² ./ num_samp
    PN_var = num_samp ./ (num_samp .-1) .* (PN²_mean - PN_mean.^2)

    return num_samp, PN_mean, PN_var, PN_hist
end
