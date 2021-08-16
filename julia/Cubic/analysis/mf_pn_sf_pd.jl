using DataFrames
using Plots
using LinearAlgebra
using StatsBase
using CSV
const modulepath =  "/Users/pcs/codes/chain/ABF2D/module"
push!(LOAD_PATH, modulepath)
using ABF2D
using Lattice
using PN
using Statistics

function readdir_csv(dir)
    files = re  addir(dir, join = true)
    filter!(s->occursin(r".csv", s), files)
end

function range_mean(x, y, r1, r2)
    idx = findall(a -> (r1 <= a < r2), x)
    return mean(y[idx])
end

## Configuration
raw_dir = "/Users/pcs/data/cubic/"
# save_dir = "/Users/pcs/data/ABF2D/gipr-pd-sf/analysis/"
sizes = 5:15:25
R = [1000 ;400; 100; 50; 15]  # Realizations
q = range(-3, 3, length = 31)
## Read all filenames
fn = Dict()
fn_temp = readdir_csv(raw_dir)
for i in 1:length(sizes)
    temp_s1=split(fn_temp[i],"L")
    temp_s2=split(temp_s1[2],"_")
    fn[parse(Int, temp_s2[1])] = fn_temp[i]
end
## Process the data
P_mean = Array{Float64}(undef, length(q),length(sizes))
P_std = similar(P_mean)
out_mean = Dict()
for i in 1:length(sizes)
    df = CSV.read(fn[sizes[i]], DataFrame)
    E_min = df.E[1]
    E_max = df.E[end]
    E_div = range(E_min, E_max, length = 20)
    out_mean[sizes[i]] = Array{Float64}(undef, length(q), R[i])
    for r in 1:R[i]
        df_r = df[df.R .== r, :]
        r1 = E_div[11]
        r2 = E_div[12]
        out_mean[sizes[i]][:,r] = [mean(df_r[(df_r.E .>= r1) .& (df_r.E .< r2), i+1]) for i in 1:length(q)]
    end
    P_mean[:, i] = mean(out_mean[sizes[i]], dims = 2)
    P_std[:, i] = std(out_mean[sizes[i]], dims = 2)
end
P_mean = collect(transpose(P_mean))
p = plot();
scatter!(p, log10.(sizes), log10.(P_mean), legend = false);
display(p)

## Linear Fitting
N = length(sizes)
X = hcat(ones(N), log10.(sizes))
Y = log10.(P_mean)
v = X\Y # v[1]: y-intercept, v[2]: slope

# se = reduce(hcat,[log10.(P_mean[:,i]) .- (v[1,i] .+ v[2,i]*log10.(sizes_int)) for i in 1:10])
# σ_y = sqrt.(1/(N-2)*vec(sum(x->x^2, se, dims = 1)))
# Δ = N*sum(x->x^2,log10.(sizes_int)) - sum(log10.(sizes_int))^2
#
# err_slope = σ_y*sqrt(N/Δ)

## Plotting
p = plot(framestyle = :box, xlabel = "L", ylabel = "GIPR", dpi = 200, grid = false);

for i in 1:length(q)
    plot!(p, log10.(sizes), v[1,i] .+ v[2,i]*log10.(sizes), label  = :none, palette = palette([:blue, :green], 31))
end
for i in 1:length(q)
    scatter!(p, log10.(sizes), log10.(P_mean[:,i]), label = string.(q), palette = palette([:blue, :green], 31))
end
xticks!(p, log10.(sizes), string.(sizes));
display(p)
##

p1 = plot(vec(q[22:31]), vec(v[2,22:31]),  c = :black,
    line = (:dot, 2), marker = (:c, 4, 0.9, Plots.stroke(0, :b)),
    label = "Our result", lw =2, frame = :box, dpi = 300, grid = false)
    plot!(p1, vec(q[22:31]), vec(v[2,22:31]), label = :none, c = :black, lw = 0, ms = 5)
    plot!(p1, vec(q[22:31]), -3*(vec(q[22:31]) .- 1),
    line = (:dot, 2), marker = (:c, 4, 0.9, Plots.stroke(0, :r)),
    c = :red, label = "Normal scaling", lw = 2)
    xlabel!("q")
    ylabel!("τ(q), Exponent of GIPR(q)")
    display(p1)


p2 = plot(vec(q), -vec(v[2,:]) .- 2*(vec(q) .- 1), c = :black,
    line = (:dot, 2), marker = (:c, 4, 0.8, Plots.stroke(0, :b)),
    legend = false, frame = :box, dpi = 300, grid = false)
    plot!(p2, vec(q), -vec(v[2,:]) .- 2*(vec(q) .- 1), yerror = mse, label = :none, c = :black, lw = 0, ms = 5)

    xlabel!("q")
    ylabel!("Anomalous exponent")
## Legendre Transformation
τ_q = -vec(v[2,:])
q = vec(q)
α = diff(τ_q1)./(q[2] - q[1])
f_α = [q[i]*α1[i] - τ_q[i] for i in 1:9]
xerr = [sqrt(mse[i+1]^2 + mse[i]^2)/(q[2] - q[1]) for i in 1:9]
yerr = err_slope

p3 = plot(α, f_α,  c = :black,
    line = (:dot, 2), marker = (:c, 4, 0.9, Plots.stroke(0, :b)),
    legend = false, frame = :box, dpi = 300, grid = false)
    plot!(p3, α, f_α, xerror = xerr, yerror = yerr, label = :none, c = :black, lw = 0, ms = 5)
    xlabel!("α")
    ylabel!("f(α), Multifractal Spectrum")

savefig(p, save_dir*"fig1.png")
savefig(p1, save_dir*"fig2.png")
savefig(p2, save_dir*"fig3.png")
savefig(p3  , save_dir*"fig4.png")


println()
println(v[2])
