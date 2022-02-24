using ABF
using Lattice
using PN
using LinearAlgebra
using SparseArrays
using Random
using Plots
using Arpack
using Statistics
using CSV, DataFrames
using LaTeXStrings
dir = "$(homedir())/2d-hop-mf-data.csv"
df = CSV.read(dir, DataFrame)
marker = (:circle, 3., 1., stroke(0, 1., :black))
line = (:line, :solid, 1.)
default(
    framestyle = :box,
    size = (600,400),
    # right_margin = [3mm 0mm],
    grid = false,
    minorticks = true,
    legend = (0.1, 0.72),
    fontfamily = "computer modern",
    tickfontsize = 14,
    guidefontsize = 16,
    legendfontsize = 14,
    annotationfontsize = 14, palette = :default)


function filt(x, y)
    return x == y
end

L = 20:20.:140.
p = plot()
for i in 1:length(L)
    dft = df[filt.(df.L, L[i]), :]
    plot!(p, dft.a, dft.fa)
end
display(p)

p = plot(legend = :topright, palette =:default)
xlabel!(L"q");
ylabel!(L"\alpha");
for i in 1:length(L)
    dft = df[filt.(df.L, L[i]), :]
    scatter!(p, dft.q, dft.a, line = line, marker = marker, label = L"$L= %$(string(convert(Int64, dft.L[1])))$")
end
display(p)

p = plot(legend = :topleft, palette =:default)
xlabel!(L"q")
ylabel!(L"\tau_q")
for i in 1:length(L)
    dft = df[filt.(df.L, L[i]), :]
    scatter!(p, dft.q, dft.tau, line = line, marker = marker, label = L"$L= %$(string(convert(Int64, dft.L[1])))$")
end
plot!(p,dft.q, 2*(dft.q .- 1), c = :black, label = "metal")
hline!([0], c = "black", line = :dot, label = "insulator")

display(p)
savefig(p, "$(homedir())/tempdata/hd.pdf")
