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
    files = readdir(dir, join = true)
    filter!(s->occursin(r".csv", s), files)
end

function cut_mean(x, y, r)
    idx = findall(a -> abs(a) < r*x[end], x)
    return mean(y[idx])
end

## Configuration
raw_dir = "/Users/pcs/data/ABF2D/gipr-pd-sf/"
save_dir = "/Users/pcs/data/ABF2D/gipr-pd-sf/analysis/"
sizes = ["L20" "L40" "L60" "L80" "L100" "L120" "L140"]
q = [1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0]
sizedir = raw_dir.*sizes
sizes_int = collect(20:20:140)
## Read all filenames
fns = Dict()
for i in 1:length(sizes)
    fns[sizes[i]] = readdir_csv(sizedir[i])
end
## Process the data
r = 0.05 # Relative energy width
P_mean = Dict() # mean of pn for each realizations
P_mean_mean = Array{Float64}(undef, length(sizes), length(q)) # ensemble average of mean of PN_mean
P_mean_std = similar(P_mean_mean)

for i in 1:length(sizes)
    P_mean[sizes_int[i]] = Float64[]
    for j in 1:length(fns[sizes[i]])
        df = CSV.read(fns[sizes[i]][j], DataFrame)
        df.E .-= 1.
        row_temp = reshape([cut_mean(df.E, df[:,i], r) for i in 2:size(df,2)], 1, size(df,2) - 1)

        if j == 1
            P_mean[sizes_int[i]] = row_temp
        else
            P_mean[sizes_int[i]] = vcat(P_mean[sizes_int[i]], row_temp)
        end
    P_mean_mean[i, :] =  mean(P_mean[sizes_int[i]], dims = 1)
    P_mean_std[i, :] = std(P_mean[sizes_int[i]], dims = 1)
    end
end


display(P_mean_mean)

## Linear Fitting
N = length(sizes)
X = hcat(ones(N), log10.(sizes_int))
Y = log10.(P_mean_mean)
v = X\Y # v[1]: y-intercept, v[2]: slope
se = reduce(hcat,[log10.(P_mean_mean[:,i]) .- (v[1,i] .+ v[2,i]*log10.(sizes_int)) for i in 1:10])
σ_y = sqrt.(1/(N-2)*vec(sum(x->x^2, se, dims = 1)))
Δ = N*sum(x->x^2,log10.(sizes_int)) - sum(log10.(sizes_int))^2

err_slope = σ_y*sqrt(N/Δ)

## Linear Fitting of exponent
X = vec(q)
Y = hcat(ones(length(q)),vec(v[2,:]))
exp_q = X\Y

## Plotting
p = plot(framestyle = :box,xaxis = :log10 , yaxis=:log10, xlabel = "L", ylabel = "GIPR", dpi = 300, grid = false)
    xticks!(20:20:140,string.(collect(20:20:140)))
    for i in 1:length(q)
        plot!(p, sizes_int, 10 .^(v[1,i] .+ v[2,i]*log10.(sizes_int)), legend = false, ls = :dot, palette = :Dark2_5, label = :none, linewidth = 2)
    end
        scatter!(p, sizes_int, P_mean_mean,
        legend = :outertopright, label = "q = ".*string.(q), yerror = P_mean_std)



p1 = plot(vec(q), vec(v[2,:]),  c = :black,
    line = (:dot, 2), marker = (:c, 4, 0.9, Plots.stroke(0, :b)),
    label = "Our result", lw =2, frame = :box, dpi = 300, grid = false)
    plot!(p1, vec(q), vec(v[2,:]), yerror = err_slope, label = :none, c = :black, lw = 0, ms = 5)
    plot!(p1, vec(q), -2*(vec(q) .- 1),
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
