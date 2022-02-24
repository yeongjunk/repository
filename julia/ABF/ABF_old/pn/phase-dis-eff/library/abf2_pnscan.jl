using Random
using LinearAlgebra
using StatsBase
include("abf2_ham.jl")
include("abf2_cores.jl")

@doc """
A struct where parameters are stored. This is the input parameters for the scanning.
    θ is in unit of π rad.
"""
struct Parameters{F<:AbstractFloat,N<:Int64, B<:Bool}
    N::N # Number of unit cell.
    W::F # If you set it to 0, only phase disorder will be considered. Otherwise, it does not affect.

    #Time reversal breaking fraction.
    V_min::F
    V_max::F
    V_num::N
    V_log::B # If true, V = 10^V

    #Manifold angle scan parameters. θ = θI = θII
    θ_min::F
    θ_max::F
    θ_num::N


    seed::N
    R::N    # Number of realizations.

    #PN bin histogram
    PN_bin_num::N

    #E bin for histogram
    E_leftedge::F
    E_rightedge::F
    E_bin_num::N
    E_abs_bw::B

    core::N # 1: SF, on+phase, 2: SF, only phase, 3: FD, on+phase, 4: FD, only phase
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

        config["PN_bin_num"],

        config["E_leftedge"],
        config["E_rightedge"],
        config["E_bin_num"],
        config["E_abs_bw"],
        config["core"]
        )

    return p
end

@doc """
Expand domain parameters
"""
function expand_params(p::Parameters)
    θ = collect(range(p.θ_min, p.θ_max, length = p.θ_num))
    V = collect(range(p.V_min, p.V_max, length = p.V_num))
    if p.V_log
        V = 10 .^V
    end
    E = collect(midpoints(range(p.E_leftedge, p.E_rightedge, length = p.E_bin_num+1)))
    return θ, V, E
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
            j += 1
        end
        for i in 2:length(x_edges)
            while (j <= length(x)) && (x_edges[i-1] <= x[j] < x_edges[i])
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
            while (j <= length(x)) && (x_edges[i-1] <= x[j] < x_edges[i])
                y_binned_sum[i-1] += y[j]
                j += 1
            end
        end
        return y_binned_sum
    end
end

#----------------------------------------------#

@doc """
It prescan 1000 times and performs the followings.
1. Scan the bandwidth of the eigen energies, which is needed to normalize the eigenvalues
2. Estimate the maximum/minimum of edges of PN bins of the histogram
"""
function prescan(f::Function, p::Parameters)
    θ, V, _ = expand_params(p)
    out_size = (p.V_num, p.θ_num)
    E_bw = zeros(Float64, out_size)

    PN_edges = Array{Float64}(undef, (out_size...,2))
    ΣPN = zeros(Float64, out_size)
    ΣPN² = zeros(Float64, out_size)
    num_samp = zeros(Int64, out_size)

    H_FD = ham_FD(2, p.N)
    T = uc_redef(p.N)
    U = LUT.(θ, p.N) # Unitaries
    H_FE = similar(U)

    for i in 1:length(θ)
        U[i] = U[i]*T*U[i] #unitary manifold
        H_FE[i] = U[i]*H_FD*U[i]'  #hamiltonain manifold
    end
    for n = 1:100 # Enough for bandwidth estimation
        rng = MersenneTwister(p.seed + n)
        r_on =  rand(rng, 2p.N) .- 0.5
        r_off = rand(rng, 12p.N) .- 0.5

        D_dis = spdiagm(0 => r_on) # onsite disorder independent of θ.
        @Threads.threads for j = 1:length(θ)
            T_dis = off_phase_dis2(H_FE[j], r_off)
            for i = 1:length(V)
                E, PN = f(H_FE[j], D_dis, T_dis, U[j], p.W, V[i])
                E_bw[i,j] = maximum(E) - minimum(E) > E_bw[i,j] ? maximum(E) - minimum(E) : E_bw[i,j]

                num_samp[i, j] += length(E)
                ΣPN[i, j] += sum(PN)
                ΣPN²[i, j] += sum(PN.^2)
            end
        end
    end

    PN_mean = ΣPN./num_samp
    PN_var = (ΣPN .- ΣPN² ./ num_samp) ./ (num_samp .- 1)
    PN_edge_hw = 4*sqrt.(PN_var)
    PN_edges[:, :, 2] .= PN_mean .+ PN_edge_hw
    for i in 1:size(PN_edges,1), j in 1:size(PN_edges,2)
        PN_edges[i, j, 1] = (PN_mean[i,j] - PN_edge_hw[i, j]) > 1. ? PN_mean[i, j] - PN_edge_hw[i, j] : 1.
    end
    return E_bw, PN_edges
end

@doc """
Check if the parameter include computing clean ABF.
"""
function isclean(p::Parameters)
    θ, V, _ = expand_params(p)
    return p.W .== 0 && p.V_log == false && any(V .== 0)
end



@doc """
Scan histogram, mean, variance of PN of R, scanning W, θ.
"""
function scan_pn(f::Function, p::Parameters)
    #-------------------Initialize output array------------------------------------#
    θ, V, _ = expand_params(p)
    out_size = (p.E_bin_num, p.V_num, p.θ_num)
    num_samp = zeros(Int64, out_size)
    ΣPN  = zeros(Float64, out_size)
    ΣPN² = zeros(Float64, out_size)
    #-------------------Prescan----------------------------------------------------#
    println("Start bandwidth scan.")
    E_bw, PN_bw = prescan(f, p)
    #-------------------Initialize Histograms--------------------------------------#
    PN_hist = Array{Histogram}(undef, p.V_num, p.θ_num)
    E_edges = range(p.E_leftedge, p.E_rightedge, length = p.E_bin_num + 1)
    PN_edges = range.(PN_bw[:,:,1], PN_bw[:,:,2], length = p.PN_bin_num + 1)
    for i in 1:p.V_num, j in 1:p.θ_num
        PN_hist[i,j] =  fit(Histogram, ([],[]), (E_edges, PN_edges[i, j]))
    end
    #-------------------Prepare Hamiltonians and Unitaries-------------------------#
    H_FD = ham_FD(2, p.N)
    T = uc_redef(p.N)
    U = LUT.(θ, p.N)
    H_FE = similar(U)
    for j in 1:length(θ)
        U[j] .= U[j]*T*U[j] #unitary manifold
        H_FE[j] = U[j]*H_FD*U[j]'  #hamiltonain manifold
    end
    #-------------------Start------------------------------------------------------#
    println("Start PN scan.")
    for n = 1:p.R
        rng = MersenneTwister(p.seed + n)
        r_on =  rand(rng, 2p.N) .- 0.5
        r_off = rand(rng, 12p.N) .- 0.5
        D_dis = spdiagm(0 => r_on) # onsite disorder independent of θ.

        @Threads.threads for j = 1:length(θ)
            T_dis = off_phase_dis2(H_FE[j], r_off)
            for i = 1:length(V)
                #-------------------core-------------------------------------------------------#
                E, PN = f(H_FE[j], D_dis, T_dis, U[j], p.W, V[i], bw = E_bw[i,j])
                #-------------------Temporary output-------------------------------------------#
                num_temp, ΣPN_temp = binned_sum(E, PN, E_edges, count = true)
                ΣPN²_temp = binned_sum(E, PN.^2, E_edges)
                PN_hist_temp = fit(Histogram, (E, PN), (E_edges, PN_edges[i,j]))
                #-------------------Merges output----------------------------------------------#
                merge!(PN_hist[i,j], PN_hist_temp)
                num_samp[:, i, j] .+= num_temp
                ΣPN[:, i, j] .+= ΣPN_temp
                ΣPN²[:, i, j] .+= ΣPN²_temp
            end
        end
    end

    PN_mean = ΣPN ./ num_samp
    PN²_mean = ΣPN² ./ num_samp
    PN_var = num_samp ./ (num_samp .-1) .* (PN²_mean .- PN_mean.^2)

    return num_samp, PN_mean, PN_var, PN_hist, E_bw
end
