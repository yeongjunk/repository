# using Plots
using LinearAlgebra
using Plots
using Lattice
using ABF
using PN
using MFA_BC
using DelimitedFiles

loaddir = "/Users/pcs/data/ABF3D/partial/partial/"

function readdir_csv(dir, join = true)
    files = readdir(dir, join = join)
    filter!(s->occursin(r".csv", s), files)
end
L = 80
q = range(2, 2, length = 1)
b = [10; 8; 5; 4;2]

fns = readdir_csv(loaddir)
idx = vec([7 3 4 5 6 13 8 9 10 11 12])
fns = fns[idx]
α = Array{Float64}[]
ltc = Lattice3D(L, L, L, 1)
bc_params = BC_Params(ltc, b, q)

for i in 1:length(fns)
    data = readdlm(fns[i], ',', Float64, '\n')
    E = data[1,:]
    psi = data[2:end, :]

    P_q = zeros(Float64, length(b), length(q))
    for j = 1:10 #αth energy
        E_r = E[j]
        ψ_r = psi[:,j]
        P_q .+= gipr(ψ_r, bc_params)
    end
    P_q /= 10

    α_temp = diff(log.(vec(P_q.^-1)))./diff(log.(L./b))
    push!(α, α_temp)
end

p = scatter()
for i in 1:length(α)
    scatter!(p, i*ones(length(α[i])), α[i], label = :none, c = 1:4, palette = :Greens_5)
end

p
