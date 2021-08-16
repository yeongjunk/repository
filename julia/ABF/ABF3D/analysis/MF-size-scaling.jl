using DataFrames
using Plots
using LinearAlgebra
using StatsBase
using JLD2
const modulepath =  "/Users/pcs/codes/chain/ABF3D/module"
push!(LOAD_PATH, modulepath)
using ABF3D
using Lattice
using PN
using Statistics

loaddir = "/Users/pcs/data/ABF3D/psi-midband"

function readdir_jld(dir, join = true)
    files = readdir(dir, join = join)
    filter!(s->occursin(r".jld", s), files)
end
L = [8 16 32]; R = [100 30 15]
q = range(-3, 3, length = 30)
b = [4; 8; 16]

fns = readdir_jld(loaddir); fns = fns[[3; 1; 2]] # sort
##
out = zeros(Float64, length(L), length(q))

for i in 1:3
    file = jldopen(fns[i], "r")
    ltc = Lattice3D(L[i],L[i],L[i],1)

    for r in 1:R[i]
        E_r = real.(file["data"][r][1])
        ψ_r = file["data"][r][2]
        for j in 1:10
            out[i, :] .+= vec(GIPR3D(ltc, ψ_r[:,j], 1, q))
        end
    end
    out[i, :] ./= (10R[i])
end
L = vec(L)
X = [ones(length(b)) log10.(L)]
Y = log10.(out)

v = X\Y

p = plot(framestyle = :box, xlabel = "L", ylabel = "GIPR", dpi = 300, grid = false);
for i in 1:length(q)
    scatter!(p, log10.(L), log10.(out[:,i]), legend = false, palette = cgrad(:viridis, 31, categorical = true, scale = log))
end
for i in 1:length(q)
    plot!(p, log10.(L), (v[1,i] .+ v[2,i]*(log10.(L))), legend = false, ls = :dot, palette = cgrad(:viridis, 31, categorical = true, scale = log), label = :none, linewidth = 2)
end

display(p)

p1 = plot(vec(q), vec(v[2,:]),  c = :black,
    line = (:dot, 2), marker = (:c, 4, 0.9, Plots.stroke(0, :b)),
    label = "Simulation", lw =2, frame = :box, dpi = 300, grid = false);
plot!(p1, vec(q), -3*vec(q .- 1), line = (:dot, 2), marker = (:c, 4, 0.9, Plots.stroke(0, :r)),
    c = :red, label = "Metallic(Theory)", lw = 2);

xlabel!("q");
ylabel!("τ(q)");
display(p1)


τ_q = -vec(v[2,:])
q = vec(q)
α = diff(τ_q)./(q[2] - q[1])
f_α = [q[i]*α[i] - τ_q[i] for i in 1:length(q)-1]

p3 = plot(α, f_α,  c = :black,
    line = (:dot, 2), marker = (:c, 4, 0.9, Plots.stroke(0, :b)),
    legend = false, frame = :box, dpi = 300, grid = false)
    xlabel!("α")
    ylabel!("f(α), Multifractal Spectrum")
