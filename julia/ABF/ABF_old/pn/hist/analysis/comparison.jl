using MAT
using Plots
using Dates

function plot_PN(W_ipr,θ_ipr, ipr, idx)
    plot(W_ipr[2:end], ipr[idx,2:end], legend = false)
    xlabel!("W")
    ylabel!("PN")
    title!("θ = $(θ_ipr[idx])"*"pi")
end

function plot_TMM(W_tmm,θ_tmm, ξ, idx)
    plot(W_tmm, ξ[:,idx], legend = false)
    xlabel!("W")
    ylabel!("ξ")
    title!("θ = $(round(θ_tmm[idx]/pi, digits = 3))"*"pi")
end


include("../library/abf2_iprscanlib.jl")
dir = "/Users/pcs/Data/ABF/part/E1"
cd(dir)
datestring =Dates.format(now(), "ddmmyy-HHmmss")

file = "/part.mat"
vars = matread(dir*file) #Dictionary

file2 = "/Users/pcs/Data/ABF/part/tmm.mat"
vars2 = matread(file2)
ipr = vars["IPR"]
ξ = vars2["xi_all"]
ξ_fb_diag = Array{Float64}(undef, 51,51)

for i in 1:51, j in 1:51
    ξ_fb_diag[j,i] = ξ[1,j,i,i]
end
W_ipr = collect(LinRange(vars["W_min"], vars["W_max"], vars["W_num"]))
W_tmm = collect(LinRange(vars2["W_min"], vars2["W_max"], vars2["W_num"]))
θ_tmm = collect(LinRange(vars2["theta_min"][1], vars2["theta_max"][1], vars2["theta_num"][1]))
θ_ipr = collect(LinRange(vars["theta_min"][1], vars["theta_max"][1], vars["theta_num"][1]))

i = 51
plot_PN(W_ipr,θ_ipr, ipr, i)
(val, idx) = findmin(abs.(θ_tmm .- θ_ipr[i]π))
plot_TMM(W_tmm,θ_tmm,ξ_fb_diag,idx)

p = config_read(vars)
print(p)
