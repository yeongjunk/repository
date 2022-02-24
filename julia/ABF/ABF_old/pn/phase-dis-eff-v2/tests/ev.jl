using Plots
using Random
using LinearAlgebra
include("../library/abf2_cores.jl")
include("../library/abf2_ham.jl")
include("../library/abf2_pnscan.jl")

#GIF Save directory

#-------------------Parameters-------------------#
θ = 0.125;
N = 3;
W = 0.00000001;
V = 1.0;

#-------------------Initialize-------------------#
r_on = (rand(2N) .- 0.5)
r_off =(rand(12N) .- 0.5)


UI = LUT(θ, N)
T = uc_redef(N)
UII = LUT(θ, N)
U = UI * T * UII

#Entangle
H_FD = ham_FD(1, N)
H_FE = U * H_FD * U'

D = spdiagm(0 => r_on)
T_dis = off_phase_dis2(H_FE, r_off)
H_SF = ham_sf_ph(H_FE, D, T_dis, U, W, V)
eig = eigen(Hermitian(Array(H_SF^2)))

n = 200
p1 = plot(real.(eig.vectors[:, n]))
plot!(p1, imag.(eig.vectors[:, n]))
plot!(p1, abs.(eig.vectors[:, n]))

p2 = plot(log.(abs.(real.(eig.vectors[:, n]))))
plot!(p2, log.(abs.(imag.(eig.vectors[:, n]))))


plot(p1, p2)
