using MAT
using LinearAlgebra
using Plots
using Statistics
using Glob
using LsqFit
dir = "/Users/pcs/data/cubic/"
ext = "*.mat"
files = glob(ext, dir)

M = 3:1:10
W = matread(files[1])["W"]
xi = [matread(files[i])["xi"] for i in 1:10]
xi = reduce((x...) -> cat(x..., dims = 3), xi)

for i in 1:length(M)
    xi[:, i, :] ./= M[i]
end
xi_mean = dropdims(mean(xi, dims = 3), dims = 3)
xi_std = dropdims(std(xi, dims = 3), dims = 3)
xi_err = xi_std./sqrt(10)

p0 = plot(W, xi_mean, yerror = xi_err)
p1 = plot(W, xi_err./xi_mean*100)

p2 = plot(W, xi_mean, yerror = xi_err)
xlims!(p2, 16, 17)
ylims!(p2, 0.55, 0.65)


## Functions for fitting
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

function model_cent(x, p)
    w = shift_norm.(x[:, 1], p[1])
    u_1 = series_1.(w, Ref([1, p[16], 0]), start_zero = false)
    u_2 = series_1.(w, Ref([1]), start_zero = true)
    ϕ_1 = @. u_1*x[:, 2]^p[2]
    ϕ_2 = @. u_2*x[:, 2]^p[3]
    F = series_2.(zip(ϕ_1, ϕ_2),  Ref([p[4], p[5], p[6], 0., 0., 0., p[10], 0., 0., 0., 0., 0.,]), 6, 2)
    return F
end

## Reshape data
M_vec = repeat(M, inner = length(W))
W_vec = repeat(W, outer = length(M))
ydata = collect(reshape(xi_mean, length(xi_mean)))
xdata = hcat(W_vec, M_vec)

## centeral part
W_min = 15.
W_max = 18.
idx = findall(x -> W_min < x < W_max, xdata[:, 1])
W_cent =  W[findall(x -> W_min < x < W_max, W)]
xdata_cent = xdata[idx, :]
ydata_cent = ydata[idx, :]
xi_mean_cent = xi_mean[findall(x -> W_min < x < W_max, W), :]

# Roughly known params
W_c_guess = 16.5
ν_guess = 1.6
Λ_c_guess = 0.6
slope_guess = -0.6

W_c_ran = [16. 17.]
ν_ran = [1.5 1.7]
Λ_c_ran = [0.5 0.7]
slope_ran = [-0.7 -0.4]

p = [W_c_guess, 1/ν_guess, -1., Λ_c_guess, slope_guess, 0., 0., 0., 0., 0.1, 0., 0.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
lb = [W_c_ran[1], 1/ν_ran[2],-Inf, Λ_c_ran[1], slope_ran[1], -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf] # Lower bound for params
ub = [W_c_ran[2], 1/ν_ran[1], 0, Λ_c_ran[2], slope_ran[2], Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf] # Upper bound for params
length(p)
# y_fit = model_rel(xdata_cent, p)
# y_fit = reshape(y_fit, length(W_cent), length(M))
# scatter!(W_cent, xi_mean_cent, ms = 3)
# plot!(W_cent,y_fit, c = :red)

fit = curve_fit(model_cent, xdata_cent, vec(ydata_cent), p, lower = lb, upper = ub)
y_fit = model_cent(xdata_cent, fit.param)
y_fit = reshape(y_fit, length(W_cent), length(M))
plt = plot(W_cent, y_fit)
for i in 1:8
    scatter!(plt, W_cent,xi_mean_cent[:, i], ms = 3, c = i)
end
display(plt)

plot(W_cent, y_fit .- xi_mean_cent)

display(fit.param)
println(1/fit.param[2])

w = shift_norm.(xdata_cent[:, 1], fit.param[1])
u_1 = series_1.(w, Ref([1, fit.param[16]]), start_zero = false)
u_2 = series_1.(w, Ref([1]), start_zero = true)
ϕ_1 = @. u_1*xdata_cent[:, 2]^fit.param[2]
ϕ_2 = @. u_2*xdata_cent[:, 2]^fit.param[3]
correction = reshape(fit.param[10]*ϕ_2, length(W_cent), length(M))
y_corrected = y_fit .- correction
scatter(W_cent, y_corrected)
plot(correction)
u_1_reshape = abs.(reshape(u_1, length(W_cent), length(M))).^(- 1/fit.param[2])

plt = plot()
for i in 1:length(M)
    scatter!(plt, W_cent,log10.(u_1_reshape[:, i]), shape = :auto)
end
plt
for i in 1:length(M)
    u_1_reshape[:, i] ./= M[i]
end

plot(log10.(1 ./ u_1_reshape), y_corrected)
scatter(log10.(1 ./ u_1_reshape), xi_mean_cent .- correction)
scatter(log10.(1 ./ u_1_reshape), xi_mean_cent .- correction)

# ## model_all
# p = fit.param
#
# lb = replace(p, 0 => -Inf) .- 0.1abs.(p)
# ub = replace(p, 0 => Inf).+ 0.1abs.(p)
#
# # lb = replace(p, 0 => -Inf)
# # ub = replace(p, 0 => Inf)
#
# fit = curve_fit(model_all, xdata, vec(ydata), p, lower = lb, upper = ub)
#
# y_fit = model_all(xdata, fit.param)
# y_fit = reshape(y_fit, length(W), length(M))
#
# plt = plot(W, y_fit)
# for i in 1:8
#     scatter!(plt, W, xi_mean[:, i], ms = 3, c = i, legend = false)
# end
# display(plt)
# 1/fit.param[2]
