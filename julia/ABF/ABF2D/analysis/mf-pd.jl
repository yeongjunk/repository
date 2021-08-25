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
dir_1 = "$(homedir())/2d-phase-mf-data-1.csv"
dir_2 = "$(homedir())/2d-phase-mf-data.csv"

df = CSV.read(dir_1, DataFrame)

marker = (:circle, 3., 1., stroke(0, 1., :black))
line = (:line, :solid, 1.0)
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

L = collect(20:20.:140.)
push!(L, 160.)
push!(L, 200.)

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
savefig(p, "$(homedir())/tempdata/pd.pdf")

p = plot(legend = :topleft, palette =:default)
xlabel!(L"q")
ylabel!(L"\tau_q")
for i in 1:length(L)
    dft = df[filt.(df.L, L[i]), :]
    scatter!(p, dft.q, dft.tau .- 2*(dft.q .- 1), line = line, marker = marker, label = L"$L= %$(string(convert(Int64, dft.L[1])))$")
end
plot!(p,dft.q, , c = :black, label = "metal")

display(p)
savefig(p, "$(homedir())/tempdata/pd-Delta.pdf")



results = Array{Float64}(undef, length(dft.q), length(L))
for i in 1:length(L)
    dft = df[filt.(df.L, L[i]), :]
    results[:, i] = dft.tau
end


log_L = log.(L)
log_tau = log.((results[21, :]))
X = [ones(length(log_L)) log_L]
Y = log_tau
v = X\Y

p = plot(legend = :topleft, palette =:default)
xlabel!(L"\ln L")
ylabel!(L"\ln \tau_2")

plot!(p, log_L, log_tau, label = L"\tau_2", marker = marker, line = line)
plot!(log_L, v[2]*log_L .+ v[1], ls = :dot, lw = 2, label = "Fitting")


savefig(p, "$(homedir())/tempdata/pd_scale.pdf")
