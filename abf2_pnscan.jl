using Random
using LinearAlgebra
using StatsBase
include("abf2_ham.jl")
include("abf2_disorder.jl")


@doc """
A struct where parameters are stored. This is the input parameters for the scanning.
    θ is in unit of π rad.
"""
struct Parameters{F<:AbstractFloat,N<:Int64, B<:Bool}
    N::N # Number of unit cell.
    W::F # This does not affect the result.

    #Time reversal breaking fraction.
    δ_min::F
    δ_max::F
    δ_num::N
    δ_log::B # If true, δ = 10^δ

    #Manifold angle scan parameters. θ = θI = θII
    θ_min::F
    θ_max::F
    θ_num::N


    seed::N
    R::N    # Number of realizations.

    #PN bin histogram
    PN_leftedge::N
    PN_rightedge::N
    PN_bin_num::N

    #E bin for histogram
    E_leftedge::F
    E_rightedge::F
    E_bin_num::N
end

@doc """
Read configuration JSON file, and construct a Parameters struct from it.
"""
function readconfig(config::Dict)
    p = Parameters(
        config["N"],
        config["W"],

        config["delta_min"],
        config["delta_max"],
        config["delta_num"],
        config["delta_log"],

        config["theta_min"],
        config["theta_max"],
        config["theta_num"],

        config["seed"],
        config["R"],

        config["PN_leftedge"],
        config["PN_rightedge"],
        config["PN_bin_num"],

        config["E_leftedge"],
        config["E_rightedge"],
        config["E_bin_num"],
        )

    return p
end

@doc """
Expand domain parameters
"""
function expand_params(p::Parameters)
    θ = collect(range(p.θ_min, p.θ_max, length = p.θ_num))
    δ = collect(range(p.δ_min, p.δ_max, length = p.δ_num))
    E = midpoints(range(p.E_leftedge, p.E_rightedge, length = p.E_bin_num+1))
    return θ, δ, E
end

@doc """
computate PN of an eigenvector.
"""
@inline function compute_pn(eigvect::AbstractVector{ComplexF64})
    PN = sum(x -> abs2(x)^2, eigvect)
    return 1 / PN
end

@doc """
Compute eigenvalues and eigenvectors, and then compute the PN of each eigenvectors.
returns the arrays (eig.values, eig.PN(eig.vectors)) of given matrix H.
"""
function eig_pn(H)
    eig = eigen(H)

    eig_num = length(eig.values)
    PN = similar(eig.values)
    for i in 1:eig_num
        PN[i] = compute_pn(eig.vectors[:,i])
    end
    return eig.values, PN
end

@doc """
for vectors (x, y), and edges x_edges of x, return the sum of 'y's of each bin.
x must be a vector sorted in increasing order. Binning range is defined as [edge[i], edge[i+1]).
The length of the returned array is length(x_edges)-1 since x_edges = x_bins + 1.
If the optional argument count = true, then it also returns number of samples of each bin.
"""
function binned_sum(x::AbstractVector{Float64},y::AbstractVector{Float64}, x_edges::AbstractVector{Float64}; count = false)
    if count
        y_binned_sum = zeros(Float64, length(x_edges)-1)
        num_samp = zeros(Int64, length(x_edges)-1)

        j = 1 #array index for x,y
        while x[j] < x_edges[1]
            # if an error is thrown here, check eigenenergy and energy ranges
            j += 1
        end
        for i in 2:length(x_edges)
            while (j <= length(x)) && x_edges[i-1] <= x[j] < x_edges[i]
                y_binned_sum[i-1] += y[j]
                num_samp[i-1] += 1
                j += 1
            end
        end
        return num_samp, y_binned_sum
    else
        y_binned_sum = zeros(Float64, length(x_edges)-1)
        j = 1 #array index for x,y
        while x[j] < x_edges[1]
            # if an error is thrown here, check eigenenergy and energy ranges
            j += 1
        end

        for i in 2:length(x_edges)
            while (j <= length(x)) && x_edges[i-1] <= x[j] < x_edges[i]
                y_binned_sum[i-1] += y[j]
                j += 1
            end
        end
        return y_binned_sum
    end
end

@doc """
Scan histogram, mean, variance of PN of R, scanning W, θ.
"""
function scan_pn_diag(p::Parameters)

    # Domain parameters
    θ, δ, _ = expand_params(p)

    if p.δ_log
        δ = 10 .^δ
        println("delta logscale is set.")
    end

    # Initialize output arrays
    num_samp = zeros(Int64, p.θ_num, p.δ_num, p.E_bin_num) # of samples for each W, θ
    ΣPN  = zeros(Float64, size(num_samp))
    ΣPN² = zeros(Float64, size(num_samp))

    # Initialize output histogram array
    PN_hist = Array{Histogram}(undef, p.θ_num, p.δ_num)
    E_edges = range(p.E_leftedge, p.E_rightedge, length = p.E_bin_num + 1)


    #PN scales with θ. We need to take this into account to make nice histogram.
    PN_max_scale = (p.PN_rightedge/p.θ_max)*θ .+ 4
    for i in 1:p.θ_num
        for j in 1:p.δ_num
            PN_edges = range(p.PN_leftedge, PN_max_scale[i], length = p.PN_bin_num+1)    #Note that the length of the edge is larger than the bin_num by 1
            PN_hist[i,j] =  fit(Histogram, ([],[]), (E_edges, PN_edges)) #setup bins for empty histogram
        end
    end

    H_FD = ham_FD(2, p.N)    # Generate the clean FE ham for all angles
    H_FE_all = Array{typeof(H_FD)}(undef, length(θ))

    T = uc_redef(p.N)
    U = LUT.(θ, p.N)

    for i in 1:length(θ)
        U[i] .= U[i]*T*U[i] #Full unitary
    end

    for i = 1:length(θ)
        H_FE_all[i] = U[i]*H_FD*U[i]'   # Fully entangled hamiltonain manifold
    end

    println("Start scanning.")
    for n = 1:p.R #Iteration for realizations
        rng = MersenneTwister(p.seed + n)

        r_on =  p.W.*(rand(rng, 2p.N) .- 0.5)
        r_on .= convert(Vector{ComplexF64}, r_on)
        r_off_norm = rand(rng, 12*p.N) .- 0.5

        D = spdiagm(0 => r_on)

        @Threads.threads for i = 1:length(θ)
            H_FE = H_FE_all[i]
            T_dis = off_phase_dis2(H_FE, r_off_norm)

            for j = 1:length(δ)
                H_dis = H_FE + D + δ[j]*T_dis
                H_dis = U[i]'*H_dis*U[i]                       # Fully detangle
                H_dis_proj = projection(H_dis)             # Projection, scale free

                (E, PN) = eig_pn(Hermitian(Array(H_dis_proj))) # Diagonalize & compute PN (array)

                num_temp, ΣPN_temp = binned_sum(E, PN, E_edges, count = true)
                ΣPN²[i, j, :] .+= binned_sum(E, PN.^2, E_edges)
                PN_edges = range(p.PN_leftedge, PN_max_scale[i], length = p.PN_bin_num + 1)
                PN_hist_temp = fit(Histogram, (E, PN), (E_edges, PN_edges))
                merge!(PN_hist[i,j], PN_hist_temp)
                num_samp[i, j, :] .+= num_temp
                ΣPN[i, j,:] .+= ΣPN_temp
            end
        end
    end

    PN_mean = ΣPN ./ num_samp
    PN²_mean = ΣPN² ./ num_samp
    PN_var = num_samp ./ (num_samp .-1) .* (PN²_mean .- PN_mean.^2)

    return num_samp, PN_mean, PN_var, PN_hist
end
