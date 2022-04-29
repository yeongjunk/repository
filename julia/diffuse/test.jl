using LinearAlgebra
using Statistics
using Plots

function create_H(L, alpha, beta, phi)
    I = 1:L
    O = (cospi.(2*alpha*I .+ phi) .+ cospi.(2*alpha*(I .-1) .+ phi) .+ ratio*cospi.(2*alpha*I .+ phi .+ beta) .+ ratio*cospi.(2*alpha*(I .- 1) .+ phi .+ beta))/4
    T = (cospi.(2*alpha*I .+ phi)/4 .- ratio*cospi.(2*alpha*I .+ phi .+ beta))/4
    return SymTridiagonal(O, T)
end

"""
Change the real basis vector to eigenbasis 
"""
function real_to_ev!(eigenstates, vec_real)
    vec_real .= eigenstates'*vec_real
end

"""
Time evolution of vector in eigenbasis
"""
function ev_evol!(e_eigenvalues, vec_ev, t)
    vec_ev .*= (e_eigenvalues.^t)
end

"""
Change from eigenbasis to real basis
"""
function ev_to_real!(eigenstates, vec_ev)
    vec_ev .=  eigenstates*vec_ev 
end

"""
Calculate the sigma square
"""
function real_uncert_sq(j ,vec_real)  
    j_mean_sq = sum(j.*abs2.(vec_real))^2
    j_sq_mean = sum(j.^2 .* abs2.(vec_real))
    return j_sq_mean - j_mean_sq
end

t = range(0, 20, length = 1000)
exp_t = 10 .^t

L = 5000
j = 1:L
R = 3
W = 0.5 
data = zeros(length(t), R) 
println("Sysem size : $(L)", "\nRealization : $(R), \ntime length : $(length(t)), \ndisorder strength : $W")

@time for r in 1:R
    H = create_H(L)
    H .= H + W*Diagonal(rand(L) .- 0.5) # 1D Anderson model
    eigenvalues, eigenstates = eigen!(H)
    e_eigenvalues = exp.(im*eigenvalues) # exp(iEt)
    vec_real = zeros(ComplexF64, L); vec_real[L÷2] = 1. # Initial delta function wave in real coordimate
    vec_ev = copy(vec_real)
    real_to_ev!(eigenstates, vec_ev)
    for i in 1:length(exp_t)
        ev_evol!(e_eigenvalues, vec_ev, exp_t[i])
        vec_real .= vec_ev
        ev_to_real!(eigenstates, vec_real) 
        data_temp = real_uncert_sq(j, vec_real)
        data[i, r] =  data_temp
    end
end
sigmasq_mean = mean(data, dims=2)
sigmasq_std = std(data, dims=2)
sigmasq_ste = sigmasq_std ./ sqrt(R)
plot(t, sqrt.(sigmasq_mean))
