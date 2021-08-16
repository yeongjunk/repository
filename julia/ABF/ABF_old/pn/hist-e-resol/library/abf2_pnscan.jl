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

    W_min::F
    W_max::F
    W_num::N
    W_log::B

    θI_min::F
    θI_max::F
    θI_num::N

    θII_min::F
    θII_max::F
    θII_num::N
    θ_diag::B

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
        config["W_min"],
        config["W_max"],
        config["W_num"],
        config["W_log"],
        config["theta_min"][1],
        config["theta_max"][1],
        config["theta_num"][1],
        config["theta_min"][2],
        config["theta_max"][2],
        config["theta_num"][2],
        config["theta_diag"],
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

function compute_pn(H::AbstractArray{Float64, 2}, E_min::Float64, E_max::Float64)
    eig = eigen(Symmetric(H), E_min, E_max)
    eig_num = length(eig.values)
    PN = Vector{Float64}(undef, eig_num)
    compute_pn_col!(PN, eig.vectors)
    return eig.values, PN
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



@doc """
Scan histogram, mean, variance of PN of R, scanning W, θ.
"""
function scan_pn_diag(p::Parameters; logscale::Bool)
    # setup parameters
    W, θI, _ = expand_params(p)

    if logscale
        W = 10 .^W
        println("logscale set...")
    end
    # Initialize output arrays
    N    = zeros(Int64, p.θI_num, p.W_num, p.E_bin_num) # of samples for each W, θ
    ΣPN  = zeros(Float64, p.θI_num, p.W_num, p.E_bin_num)
    ΣPN² = zeros(Float64, p.θI_num, p.W_num, p.E_bin_num)

    # initialize histogram array
    PN_hist = Array{Histogram}(undef, p.θI_num, p.W_num)
    PN_bins = range(p.PN_min, p.PN_max, length = p.PN_bin_num) # Histogram bins
    E_bins = range(-0.5, 0.5, length = p.E_bin_num)

    for i in eachindex(PN_hist)
        PN_hist[i] =  fit(Histogram, ([],[]), (E_bins, PN_bins))
    end

    # allocate space for hamiltonain
    td_num = Threads.nthreads()
    H = zeros(Float64, 2p.N, 2p.N, td_num) # individual hamiltonian matrix storage for each threads.

    for n = 1:p.R #Iteration for each disorder realizations
        rng = MersenneTwister(p.seed + n)
        r_arr = rand(rng, 2p.N) .- 0.5
        for i = 1:length(θI)
            overwrite_ham_abf2!(view(H, :, :, 1), θI[i], θI[i], p.N)
            for t = 2:td_num
                view(H, :, :, t) .= view(H, :, :, 1)
            end
            Threads.@threads for j = 1:length(W) # Iteration for W

                r_onsite = W[j]*r_arr
                H_dis = view(H, :, :, Threads.threadid())

                add_diag!(H_dis, r_onsite)
                (E, PN) = compute_pn(H_dis, -2., 0.)
                E = (E .+ 1) ./ W[j]

                PN_hist_temp = fit(Histogram, (E, PN), (E_bins, PN_bins))
                N_temp, ΣPN_temp = sum_bin(E, PN, -0.5, 0.5, p.E_bin_num)
                _, ΣPN²_temp = sum_bin(E, PN.^2, -0.5, 0.5, p.E_bin_num)

                N[i, j, :] .+= N_temp
                ΣPN[i, j,:] .+= ΣPN_temp
                ΣPN²[i, j, :] .+= ΣPN²_temp

                merge!(PN_hist[i,j], PN_hist_temp)
                sub_diag!(H_dis, r_onsite) #remove disorder for recycling
            end
        end
    end
    PN_mean = ΣPN ./ N
    PN²_mean = ΣPN² ./ N
    PN_var = N ./ (N .-1) .* (PN²_mean - PN_mean.^2) #

    return N, PN_mean, PN_var, PN_hist
end

@doc """
Scan histogram, mean, variance of PN of R, scanning W, θ.
"""
function scan_pn(p::Parameters)
    # setup parameters
    W, θI, θII = expand_params(p)

    if p.W_log
        W = 10 .^W
    end
    # Initialize output arrays
    N    = zeros(Int64, p.θI_num, p.θII_num, p.W_num, p.E_binnum) # of samples for each W, θ
    ΣPN  = zeros(Float64, p.θI_num, p.θII_num,  p.W_num, p.E_binnum)
    ΣPN² = zeros(Float64, p.θI_num, p.θII_num, p.W_num, p.E_binnum)

    # initialize histogram array
    PN_hist = Array{Histogram}(undef, p.θI_num, p.θII_num,  p.W_num)
    PN_bins = range(p.PN_min, p.PN_max, length = p.PN_binnum) # Histogram bins
    E_bins = range(-1, 1, length = p.E_binnum)

    for i in eachindex(PN_hist)
        PN_hist[i] =  fit(Histogram, ([],[]), (E_bins, PN_bins))
    end

    # allocate space for hamiltonain
    td_num = Threads.nthreads()
    H = zeros(Float64, 2p.N, 2p.N, td_num) # individual hamiltonian matrix storage for each threads.

    for n = 1:p.R #Iteration for each disorder realizations
        rng = MersenneTwister(p.seed + n)
        r_arr = rand(rng, 2p.N) .- 0.5
        for i = 1:length(θI), j = 1:length(θII)
            overwrite_ham_abf2!(view(H, :, :, t), θI[i], θII[j], p.N)
            for t = 2:td_num
                view(H,:,:,t) .= view(H,:,:,1)
            end
            Threads.@threads for k = 1:length(W) # Iteration for W

                r_onsite = W[j]*r_arr
                H_dis = view(H, :, :, Threads.threadid())

                add_diag!(H_dis, r_onsite)
                (E, PN) = compute_pn(H_dis, -2., 0.)

                PN_hist_temp = fit(Histogram, (E, PN), (E_bins, PN_bins))

                N_temp, ΣPN_temp = sum_bin(E, PN, -2., 0., p.E_binnum)
                _, ΣPN²_temp = sum_bin(E, PN.^2, -2., 0., p.E_binnum)

                N[i, j, :] .+= N_temp
                ΣPN[i, j,:] .+= ΣPN_temp
                ΣPN²[i, j, :] += ΣPN²_temp

                merge!(PN_hist[i,j], PN_hist_temp)
                sub_diag!(H_dis, r_onsite) #remove disorder for recycling

            end
        end
    end
    PN_mean = ΣPN ./ N
    PN²_mean = ΣPN² ./ N
    PN_var = N ./ (N .-1) .* (PN²_mean - PN_mean.^2) #

    return N, PN_mean, PN_var, PN_hist
end
