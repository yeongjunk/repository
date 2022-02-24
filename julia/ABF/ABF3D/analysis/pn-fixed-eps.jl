using DataFrames, CSV, MAT
using ROAG, Binning
using Statistics
using StatsBase
using Glob
using Plots
using LsqFit
using LaTeXStrings


##
L = [20 30 40 50 60]
rdir = "/Users/pcs/data/ABF-sum/raw-data/3d-sf-on-pn-roag/par-ipr-fix-eps/"
dir = [rdir*"L$(L[i])/" for i in 1:length(L)]

# k = 5
# for i in 2:7
#     init = (i-1)*3
#     for j in 1:3
#         mv(dir[k]*"$(i)/L$(L[k])_Th$(j).csv", dir[k]*"$(i)/L$(L[k])_Th$(init + j).csv")
#     end
# end


θ = range(0.08, 0.12, length = 21)

savedir = "/Users/pcs/data/ABF3D/analyzed-data/"

ipr_mean = Array{Float64}(undef, length(θ), length(L))
ipr_std = similar(ipr_mean)
ipr_ste = similar(ipr_mean)
for i in 1:length(L)
    for j in 1:length(θ)
        df = CSV.read(dir[i]*"L$(L[i])_Th$j.csv", DataFrame)
        ipr_mean[j, i] = mean(df.q2)
        ipr_std[j, i] = std(df.q2)
        ipr_ste[j, i] = ipr_std[j, i] / sqrt(length(df.E))
    end
end

τ = similar(ipr_mean)
τ_err = similar(ipr_mean)
for i in 1:length(L)
    τ[:, i] = log.(ipr_mean[:, i]) ./ log(0.1)
    τ_err[:, i] = abs.(ipr_ste[:, i] ./ ipr_mean[:,i] ./ log(0.1))
end

plot(θ, τ, yerror = τ_err, label = "L = ".*string.(L), legend = :bottomright)

## model functions
shift_norm(x, x_c) = (x - x_c)/x_c

function series_1(x, a; start_zero = false)
    u = 0.
    if start_zero
        for j in 0:length(a)-1
            u += a[j + 1]*x^j
        end
    else
        for j in 1:length(a)
            u += a[j]*x^j
        end
    end
    return u
end


function series_2(r::Tuple{Real, Real}, a::Vector{Float64},  n, m)
    u = 0.
    A = reshape(a, n, m)
    for i in 0:size(A, 1)-1, j in 0:size(A, 2)-1
        u += (A[i + 1, j + 1]*getindex(r, 1)^i*getindex(r, 2)^j)
    end
    return u
end

function model(x, p)
    th = shift_norm.(x[:, 1], p[1])
    u_1 = series_1.(th, Ref([1, p[11], 0.]), start_zero = false)
    u_2 = series_1.(th, Ref([1, 0., 0.]), start_zero = true)
    ϕ_1 = @. u_1*x[:, 2]^p[2]
    ϕ_2 = @. u_2*x[:, 2]^p[3]
    F = series_2.(zip(ϕ_1, ϕ_2),  Ref([p[4], p[5], p[6], p[7], p[8], p[9], 0., 0., p[10], 0., 0., 0.,0., 0., 0., 0.]), 8, 2)
    return F
end


## Reshape data
L_vec = repeat(vec(L), inner = length(θ))
θ_vec = repeat(θ, outer = length(L))
ydata = collect(reshape(τ, length(τ)))
ywt = 1 ./ collect(reshape(τ_err, length(τ_err))).^2
xdata = hcat(θ_vec, L_vec)
## Parameters
θ_c_guess = 0.1
ν_guess = 1.6
τ_c_guess = 1.5
slope_guess = 0.5

θ_c_ran = [0.095 0.105]
ν_ran = [1.5 1.7]
τ_c_ran = [1. 2.]
slope_ran = [0. 10.]
p = [θ_c_guess, 1/ν_guess, -2., τ_c_guess, slope_guess, 0., 0., 0., 0., 0., 0., 0.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
lb = [θ_c_ran[1], 1/ν_ran[2], -5., τ_c_ran[1], slope_ran[1], -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf] # Lower bound for params
ub = [θ_c_ran[2], 1/ν_ran[1], 0, τ_c_ran[2], slope_ran[2], Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf] # Upper bound for params
p = p[1:11]
lb = lb[1:11]
ub = ub[1:11]

##
length(p)
y_fit = model(xdata, p)
y_fit = reshape(y_fit, length(θ), length(L))
scatter(θ, τ, ms = 3)
plot!(θ, y_fit, c = :red)

fit = curve_fit(model, xdata, vec(ydata), ywt, p, lower = lb, upper = ub)
se = margin_error(fit)

# fit = curve_fit(model, xdata, vec(ydata), ywt, fit.param, lower = lb, upper = ub)
# se = margin_error(fit)
marker = (:circle, 4, 1., stroke(-0.5, 1., :black))
line = (:line, :solid, 2)
palette = :Dark2_5

default(
    framestyle = :box,
    # right_margin = [3mm 0mm],
    grid = false,
    minorticks = true,
    size = (600, 400),
    legend = (0.1, 0.75),
    fontfamily = "computer modern",
    tickfontsize = 10,
    guidefontsize = 10,
    legendfontsize = 10,
    palette = :Dark2_5
    )


y_fit = model(xdata, fit.param)
y_fit = reshape(y_fit, length(θ), length(L))

θ_norm = shift_norm.(xdata[:, 1], fit.param[1])
u_1 = series_1.(θ_norm, Ref([1, fit.param[11], 0.]), start_zero = false)
u_2 = series_1.(θ_norm, Ref([1]), start_zero = true)
ϕ_1 = @. u_1*xdata[:, 2]^fit.param[2]
ϕ_2 = @. u_2*xdata[:, 2]^fit.param[3]
correction = reshape(fit.param[10]*ϕ_2, length(θ), length(L))
y_corrected = y_fit .- correction
ydata_corrected = τ .- correction

u_1_reshape = reshape(u_1, length(θ), length(L))
xi = abs.(u_1_reshape).^(-1. /fit.param[2])


plt = plot(legend = :bottomright)
for i in 1:length(L)
    if i == 1
        plot!(plt, θ, y_corrected[:, i], c = :black, line = :dot, label = "Fitting")
    else
        plot!(plt, θ, y_corrected[:, i], c = :black, line = :dot, label = :none)
    end
    scatter!(plt, θ, ydata_corrected[:, i], marker = marker, ms = 2, label = "L = $(L[i])")
    plot!(plt, θ, ydata_corrected[:, i], label =:none, c = :transparent, yerror = τ_err)

end
ylabel!(L"\tilde{\tau}_{corrected}")
xlabel!(L"\theta/\pi")
display(plt)



plt_1 = plot(legend = :bottomright)
for i in 1:length(L)
    if i == 1
        plot!(plt_1, θ, y_fit[:, i], c = :black, line = :dot, label = "fitting")
    else
        plot!(plt_1, θ, y_fit[:, i], c = :black, line = :dot, label = :none)
    end
    scatter!(plt_1, θ, τ[:, i], marker = marker, label = "L = $(L[i])")
    plot!(plt_1, θ, τ[:, i], label =:none, c = :transparent, yerror = τ_err)
end
ylabel!(L"\tilde{\tau}")
xlabel!(L"\theta/\pi")
display(plt_1)


xi_L = similar(xi)
for i in 1:length(L)
    xi_L[:, i] = xi[:, i]./ L[i]
end

plt_2 = plot(inset = (1, bbox(0.04, 0.1, 0.35, 0.35, :bottom, :right)), legend = :topright)
scatter!(plt_2, sp = 1, log10.(xi_L), y_corrected, label ="L = ".*string.(L), ,marker = marker, palette = :default)
ylabel!(sp = 1, L"\tilde{\tau}_{corrected}")
xlabel!(sp = 1, L"\log_{10} (\xi/L)")

plot!(plt_2, sp = 2, θ, log10.(xi)[:,2], legend = false)
xlabel!(plt_2, sp = 2, L"\theta/\pi")
ylabel!(plt_2, sp = 2, L"\log_{10} \xi")


savefig(plt, savedir*"tau.pdf")
savefig(plt_1, savedir*"tau_cor.pdf")
savefig(plt_2, savedir*"xi.pdf")
savefig(plt_3, savedir*"scaling_function.pdf")
