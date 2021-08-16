# using Plots
using LinearAlgebra
using JLD2

const modulepath =  "/Users/pcs/codes/chain/ABF3D/module"
push!(LOAD_PATH, modulepath)
using Lattice
using ABF3D
using PN
using MFA_BC

loaddir = "/Users/pcs/data/ABF3D/psi-E1"

function readdir_jld(dir, join = true)
    files = readdir(dir, join = join)
    filter!(s->occursin(r".jld", s), files)
end
L = [8 16 32]; R = [100 30 15]
q = range(-3, 3, length = 30)
b = [8; 16; 32;]

fns = readdir_jld(loaddir); fns = fns[[3; 1; 2]]

i = 3

ltc = Lattice3D(L[i],L[i],L[i],1)
file = jldopen(fns[i], "r")

P_q = zeros(Float64, length(b), length(q))
for r in 1:10 # rth realizations
    α = 1 #αth energy
    E_r = real.(file["data"][r][1][α])
    ψ_r = file["data"][r][2][:,α]
    bc_params = BC_Params(ltc, b, q)
    P_q .= gipr(ψ_r, bc_params)
end
τ_q = -fit_exp_gipr(P_q, bc_params; y_intercept = true)
α, f_α = mf_spectrum(τ_q, bc_params)

p = plot(framestyle = :box, xlabel = "L", ylabel = "GIPR", dpi = 300, grid = false);
for j in 1:length(q)
    scatter!(p, log10.(L[i]./b), log10.(P_q[:,j]), legend = false, palette = cgrad(:viridis, 31, categorical = true, scale = log))
end
for j in 1:length(q)
    plot!(p, log10.(L[i]./b), -(τ_q[1,j] .+ τ_q[2,j]*(log10.(L[i]./b))), legend = false, ls = :dot, palette = cgrad(:viridis, 31, categorical = true, scale = log), label = :none, linewidth = 2)
end

display(p)

p1 = plot(q, τ_q[2,:],  c = :black, legend = :bottomright,
    line = (:dot, 2), marker = (:c, 4, 0.9, Plots.stroke(0, :b)),
    label = "Simulation", lw =2, frame = :box, dpi = 300, grid = false);
plot!(p1, q, 3*(q .-1), line = (:dot, 2), marker = (:c, 4, 0.9, Plots.stroke(0, :r)),
    c = :red, label = "Metallic(Theory)", lw = 2);

xlabel!("q");
ylabel!("τ(q)");
display(p1)


p3 = plot(α, f_α,  c = :black,
    line = (:dot, 2), marker = (:c, 4, 0.9, Plots.stroke(0, :b)),
    legend = false, frame = :box, dpi = 300, grid = false)
    xlabel!("α")
    ylabel!("f(α), Multifractal Spectrum")
