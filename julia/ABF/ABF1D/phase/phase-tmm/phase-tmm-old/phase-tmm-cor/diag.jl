using Plots
using Random
using LinearAlgebra
using StaticArrays
using ABF, Lattice, PN, SparseArrays
using ProgressBars
using ProfileView

@doc """
offdiagonal phase disorder ver. 2
"""
function hoppings(θ, p)
    c = cospi(θ)
    s = sinpi(θ)
    s4 = s^4; c4 = c^4; s2c2 = s^2*c^2
    t_1 = 2im*(p[1, 1]s4 - p[1, 1]c4 + p[2, 1]*s2c2 + p[2, 2]*s2c2 
        + p[3, 1]*c4 - p[3, 2]*s2c2 - p[4, 1]*s2c2 + p[4, 2]*s4
        - p[5, 1]*c4 - p[5, 2]*s4)*s2c2
    t_2 = 2im*(-p[2, 1] + p[3, 1] + p[4, 1] - p[5, 1])*s4*c4
    return t_1, t_2
end

@doc """
transfer matrix of nth supercell... Note that tt is hoppings in (j-1)th row, t is hoppings in (j)th row.
"""
function tm(E, t_1, t_2, t_m1, t_m2, tt_1, tt_2, tt_m1, tt_m2, vartype)
    TM = Array{Complex{vartype}}(undef, 4, 4)
    
    H_n_0 = [-E  tt_1; 
             t_m1 -E]
    
    T_n_1 = [tt_2   0.;
             t_1     t_2]
    
    T_n_m1 = [tt_m2 tt_m1; 
              0.   t_m2]
    
    TM[1:2, 1:2] = -inv(T_n_1)*H_n_0
    TM[1:2, 3:4] = -inv(T_n_1)*T_n_m1
    TM[3:4, 1:2] = I(2)
    TM[3:4, 3:4] .= 0
    return TM
end

function eig_corr(x, i, len_cor)
    y = zeros(Float64, len_cor)
    for r in 0:len_cor - 1
        y[r + 1] += abs(x[i] - x[i + r])
    end
    return y
end

"""
Correlation averaged over j and realization
"""
function cor_tmm(;E = 0., θ = 0.25, δ = 0., N = 10000, q = 4, R = 10, seed = 1234)
    vartype = BigFloat
    ltc = Lattice1D(N, 2)
    
    cutoff_wf_1 = 1
    cutoff_wf_2 = N
    len_cwf = cutoff_wf_2 - cutoff_wf_1 + 1 
    len_r = len_cwf÷2
    g_r = zeros(vartype, len_r)
    g_r_sq = zeros(vartype, len_r)
    
    for r in tqdm(1:R)
        data =  loc_length2(θ = θ, E = E, N = N, q = q, seed = seed)
        logpsi = Float64.(log.(abs.(data)))
        logpsi_sq = logpsi.^2
        for i in 1:len_r
            g_r .+= eig_corr(logpsi, i, len_r)
            g_r_sq .+= eig_corr(logpsi_sq, i, len_r)
        end
    end
    
    return g_r/R/len_r, g_r_sq/R/len_r
end


