using LinearAlgebra
using Plots
using Random

hopping(α, n) = cospi(2*α*n)

function tm(E, t_m, t_m1)
   return [0. -1.; t_m1/t_m -E/t_m]
end

function ham(N, α)
   hoppings = [hopping(α, i) for i in 1:N-1]
   return SymTridiagonal(zeros(Float64, N), hoppings)
end 
α = (sqrt(5) - 1)/2
N = 10001

H = ham(N, α)
E, vecs = eigen(H)
exp_E = exp.(im*E)
psi = Float64[] 
psi_m = vecs[1:2, end÷2] 
push!(psi, psi_m[1])
@time for i in 2:N-1
    t_1 = H[i, i+1]
    t_2 = H[i, i-1] 
    T = tm(E[end÷2], t_1, t_2)
    psi_m=T*psi_m
    push!(psi, psi_m[1])
end


plot(abs.(vecs[1:end-1,end÷2]).-abs.(psi))
