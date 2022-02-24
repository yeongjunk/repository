
#= collections of slightly modified or generalized test functions for experiments =#

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

@doc """"
Same as xi_abf2, but additionally return a bool isRinf. For testing whether R -> inf during the iteration.
"""
function ξ_abf2!(E::Float64, W::Float64, θI::Float64, θII::Float64, r_all::Array{Float64, 2},isRinf::Bool)::Float64
    N::UInt64 = size(r_all,2)                                  # number of iterations
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

    for j = 1:N                                                                 # loop over the sites/iterations
        e_aj = W*r_all[1, j]                                                    # random onsite energy
        e_bj = W*r_all[2, j]                                                    # random onsite energy
        e_pj = cos2θI + e_aj*cosθII² + e_bj*sinθII²                             # random onsite energies in
        e_fj = -cos2θI + e_bj*cosθII² + e_aj*sinθII²                            # the semi-detangled basis
        tj = t*abs(e_aj - e_bj)                                                 # the disordered hopping in the semi-detangled basis

        # R_2j-1
        R = (E - e_fj)/tj - tI/tj/R
        ΠR_k *= R
        isRinf=isinf(R)
        # R_2j
        R = (E - e_pj)/tI - tj/tI/R
        ΠR_k *= R
        isRinf=isinf(R)
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
@doc """"
Same as xi_abf2, but optimization for log is not done. Also return isRinf additionally.
"""
function ξ_abf2_unopt!(E::Float64, W::Float64, θI::Float64, θII::Float64, r_all::Array{Float64, 2},isRinf::Bool)::Float64
    N::UInt64 = size(r_all,2)                                  # number of iterations
    tI = abs(sin(2θI))
    t  = abs(sin(2θII)/2)
    R = 1.
    ΣlogR = 0.

    cos2θI = cos(2θI)
    cosθII² = cos(θII)^2
    sinθII² = sin(θII)^2

    for j = 1:N                                                         # loop over the sites/iterations
        e_aj = W*r_all[1, j]                                       # random onsite energy
        e_bj = W*r_all[2, j]                                       # random onsite energy
        e_pj = cos2θI + e_aj*cosθII² + e_bj*sinθII²                     # random onsite energies in
        e_fj = -cos2θI + e_bj*cosθII² + e_aj*sinθII²                    # the semi-detangled basis
        tj = t*abs(e_aj - e_bj)                                         # the disordered hopping in the semi-detangled basis

        # R_2j-1
        R = (E - e_fj)/tj - tI/tj/R
        ΣlogR += log(abs(R))
        isRinf=isinf(R)
        # R_2j
        R = (E - e_pj)/tI - tj/tI/R
        ΣlogR += log(abs(R))
        isRinf=isinf(R)
        #= log takes huge proportion of calculation time. We will use multiplication and log is used
        #every k_max times. The drawback of this method is that if R is too large(small xi)
        #overflow error will occur(this happens near pi/2 and 0).  =#

    end

    return N/ΣlogR
end



@doc """
Same as ξ_scan, but the calculator is now argument, and returns isRinf. designed for test purpose
"""
function ξ_scan(p::Parameters, calculator::Function; single::Bool=false)
    E_all, W_all, θI_all, θII_all, r_all = params_expand(p)
    #ξ of arrays of parameters will be calculated. This will be done parallel in multiple threads.
    isRinf = false
    if single
        ξ_all = ξ_scan_core_s!(isRinf, E_all, W_all, θI_all, θII_all, r_all, calculator::Function)
    else
        ξ_all = ξ_scan_core!(isRinf, E_all, W_all, θI_all, θII_all, r_all, calculator::Function)
    end
    return ξ_all,isRinf
end
function ξ_scan_core!(isRinf::Bool, E_all::Array{Float64}, W_all::Array{Float64}, θI_all::Array{Float64}, θII_all::Array{Float64}, r_all::Array{Float64}, calculator::Function)
    ξ_all = Array{Float64}(undef, length(E_all), length(W_all), length(θI_all), length(θII_all))
    @time Threads.@threads for i in 1:length(E_all)
        for j in 1:length(W_all), k in 1:length(θI_all), l in 1:length(θII_all)
            ξ_all[i, j, k, l] = calculator(E_all[i], W_all[j], θI_all[k], θII_all[l], r_all, isRinf)
        end
        println("$(i)")
    end
    return ξ_all
end
function ξ_scan_core_s!(isRinf::Bool, E_all::Array{Float64}, W_all::Array{Float64}, θI_all::Array{Float64}, θII_all::Array{Float64}, r_all::Array{Float64}, calculator::Function)
    ξ_all = Array{Float64}(undef, length(E_all), length(W_all), length(θI_all), length(θII_all))
    @time for i in 1:length(E_all)
        for j in 1:length(W_all), k in 1:length(θI_all), l in 1:length(θII_all)
            ξ_all[i, j, k, l] = calculator(E_all[i], W_all[j], θI_all[k], θII_all[l], r_all, isRinf)
        end
        println("$(i)")
    end
    return ξ_all
end
@doc """"
Same as xi_abf2, but instead of returning ξ, this will return an array of R_2n. Not implemented
"""
function ξ_abf2_R!(E::Float64, W::Float64, θI::Float64, θII::Float64, r_all::Array{Float64, 2}, isRinf::Bool)::Array{Float64}
    N::UInt64 = size(r_all,2)                                  # number of iterations
    tI = abs(sin(2θI))
    t  = abs(sin(2θII)/2)
    R = 1.
    ΣlogR = 0.

    cos2θI = cos(2θI)
    cosθII² = cos(θII)^2
    sinθII² = sin(θII)^2


    for j = 1:N                                                         # loop over the sites/iterations
        e_aj = W*r_all[1, j]                                       # random onsite energy
        e_bj = W*r_all[2, j]                                       # random onsite energy
        e_pj = cos2θI + e_aj*cosθII² + e_bj*sinθII²                     # random onsite energies in
        e_fj = -cos2θI + e_bj*cosθII² + e_aj*sinθII²                    # the semi-detangled basis
        tj = t*abs(e_aj - e_bj)                                         # the disordered hopping in the semi-detangled basis

        # R_2j-1
        R = (E - e_fj)/tj - tI/tj/R
        ΣlogR += log(abs(R))
        isRinf=isinf(R)
        # R_2j
        R = (E - e_pj)/tI - tj/tI/R
        ΣlogR += log(abs(R))
        isRinf=isinf(R)
        #= log takes huge proportion of calculation time. We will use multiplication and log is used
        #every k_max times. The drawback of this method is that if R is too large(small xi)
        #overflow error will occur(this happens near pi/2 and 0).  =#

    end

    return N/ΣlogR
end
@doc """"
Same as xi_abf2, but this will calculate ξ from p_n+1/p_n, not f_n+1/f_n.(Not implemented yet)
"""
function ξ_abf2_p!(E::Float64, W::Float64, θI::Float64, θII::Float64, r_all::Array{Float64, 2}, isRinf::Bool)::Array{Float64}
    N::UInt64 = size(r_all,2)                                  # number of iterations
    tI = abs(sin(2θI))
    t  = abs(sin(2θII)/2)
    R = 1.
    ΣlogR = 0.

    cos2θI = cos(2θI)
    cosθII² = cos(θII)^2
    sinθII² = sin(θII)^2


    for j = 1:N                                                         # loop over the sites/iterations
        e_aj = W*r_all[1, j]                                       # random onsite energy
        e_bj = W*r_all[2, j]                                       # random onsite energy
        e_pj = cos2θI + e_aj*cosθII² + e_bj*sinθII²                     # random onsite energies in
        e_fj = -cos2θI + e_bj*cosθII² + e_aj*sinθII²                    # the semi-detangled basis
        tj = t*abs(e_aj - e_bj)                                         # the disordered hopping in the semi-detangled basis

        # R_2j-1
        R = (E - e_fj)/tj - tI/tj/R
        ΣlogR += log(abs(R))
        isRinf=isinf(R)
        # R_2j
        R = (E - e_pj)/tI - tj/tI/R
        ΣlogR += log(abs(R))
        isRinf=isinf(R)
        #= log takes huge proportion of calculation time. We will use multiplication and log is used
        #every k_max times. The drawback of this method is that if R is too large(small xi)
        #overflow error will occur(this happens near pi/2 and 0).  =#

    end

    return N/ΣlogR
end
