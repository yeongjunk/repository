using MAT
using Random
using Distributed


@doc """
A struct where parameters are stored. This will be used as inputs.
"""
struct Parameters
    E_min::Float64
    E_max::Float64
    E_num::Int64

    W_min::Float64
    W_max::Float64
    W_num::Int64

    θI_min::Float64
    θI_max::Float64
    θI_num::Int64

    θII_min::Float64
    θII_max::Float64
    θII_num::Int64

    seed::Int64
    N::Int64
end

@doc """
Read configuration JSON file, and construct a Parameters struct from it.
"""
function config_read(config::Dict{String,Any})
    p = Parameters(config["E_min"],
           config["E_max"],
           config["E_num"],

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
           config["N"],
           )

           return p
end

function params_expand(p::Parameters)
    E_all   = collect(LinRange(p.E_min, p.E_max, p.E_num))     #   Range of eigenenergies
    W_all   = collect(LinRange(p.W_min, p.W_max, p.W_num))     #   Range of disorder strengths

    #   Ranges of angles parameterising the nu=2 ABF models
    θI_all  = collect(LinRange(p.θI_min, p.θI_max, p.θI_num));
    θII_all  = collect(LinRange(p.θII_min, p.θII_max, p.θI_num));

    r_all   = rand_array_create(p.seed, p.N)              # pregenerate random numbers

    return E_all, W_all, θI_all, θII_all, r_all
end

function params_iter(p::Parameters)
    E_all   = LinRange(p.E_min, p.E_max, p.E_num)       #   Range of eigenenergies
    W_all   = LinRange(p.W_min, p.W_max, p.W_num)     #   Range of disorder strengths

    #   Ranges of angles parameterising the nu=2 ABF models
    θI_all  = LinRange(p.θI_min, p.θI_max, p.θI_num);
    θII_all  = LinRange(p.θII_min, p.θII_max, p.θI_num);

    r_all   = rand_array_create(p.seed, p.N)              # pregenerate random numbers
    return E_all, W_all, θI_all,θII_all, r_all
end

@doc """
read ξ and parameters from saved .mat file
"""
function ξ_read(matfilename::String)
    vars = matread(matfilename) #Dictionary

    # Localization length
    ξ_all   = vars["xi_all"]
    params = config_read(vars) #create a Parameters struct

    return ξ_all, params
end


@doc """
Pregenerate random numbers needed for the computation of the localisation
length
"""
function rand_array_create(seed::Int64, N::Int64)
    rng = MersenneTwister(seed)                                         # initialise RNG
    r_all = Array{Float64}(undef, (2, N))                                #
    rand!(rng, r_all)                                          # populate the array by random numbers in [0, 1]
    r_all .-= 0.5
    # shift the values to [0.5, 0.5]

    return r_all
end

@doc """
Compute localisation length for nu=2 ABF model
"""
function ξ_abf2(E::Float64, W::Float64, θI::Float64, θII::Float64, r_all::Array{Float64, 2})::Float64
    N::UInt64 = size(r_all,2)                                   # number of iterations

    tI = abs(sin(2θI))
    t  = abs(sin(2θII)/2)
    R = 1.
    ΣlogR = 0.
    ΠR_k = 1.

    cos2θI = cos(2θI)
    cosθII² = cos(θII)^2
    sinθII² = sin(θII)^2

    k = 1
    k_max = 10

    for j = 1:N                                                         # loop over the sites/iterations
        e_aj = W*r_all[1, j]                                       # random onsite energy
        e_bj = W*r_all[2, j]                                       # random onsite energy
        e_pj = cos2θI + e_aj*cosθII² + e_bj*sinθII²                     # random onsite energies in
        e_fj = -cos2θI + e_bj*cosθII² + e_aj*sinθII²                    # the semi-detangled basis
        tj = t*abs(e_aj - e_bj)                                         # the disordered hopping in the semi-detangled basis

        # R_2j-1
        R = (E - e_fj)/tj - tI/tj/R
        ΠR_k *= R

        # R_2j
        R = (E - e_pj)/tI - tj/tI/R
        ΠR_k *= R

        #= log takes huge proportion of calculation time. We will use multiplication and log is used
        #every k_max times. The drawback of this method is that if R is too large(small xi)
        #overflow error will occur(this happens near pi/2 and 0).  =#

        if k == k_max
            ΣlogR += log(abs(ΠR_k))
            ΠR_k = 1.
            k = 1
        else
            k += 1
        end

    end

    return N/ΣlogR
end

@doc """
(Distributed)Compute localisation length by scanning the parameters. Version 1(using distributed)
"""
function ξ_scan_dist(p::Parameters)
    E_all, W_all, θI_all, θII_all, r_all = params_expand(p)
    return ξ_scan_core_dist(E_all, W_all, θI_all, θII_all, r_all)
end

function ξ_scan_core_dist(E_all::Array{Float64}, W_all::Array{Float64}, θI_all::Array{Float64}, θII_all::Array{Float64}, r_all::Array{Float64})
    ξ_all = Array{Float64}(undef, length(E_all), length(W_all), length(θI_all), length(θII_all))
    len_div = length(θII_all)
    #By combining nested for loop into one for loop distributed computing with # cores > 100 is possible
    @time @sync @distributed for k in 1:length(θII_all)*length(θI_all)
        @fetchfrom 1 len_div
        x = div(k-1,len_div)+1
        y = mod1(k,len_div)
        for i in 1:length(E_all), j in 1:length(W_all)
            ξ_all[i, j, x, y] = ξ_abf2(E_all[i], W_all[j], θI_all[x], θII_all[y], r_all)
        end
    end
    return ξ_all
end
@doc """
(Distributed)Compute localisation length by scanning the parameters. Version 2(using pmap)
"""

function ξ_scan_pmap(p::Parameters)
    ξ_all = Array{Float64}(undef, p.E_num, p.W_num, p.θI_num, p.θII_num)
    E_all, W_all, θI_all, θII_all, r_all = params_iter(p)
    ξ_all[:,:,:,:] = pmap(t->((E,W,θI,θII)=t; ξ_abf2(E, W, θI, θII, r_all)), Iterators.product(E_all, W_all, θI_all,θII_all), batch_size = 50)
    return ξ_all
end
