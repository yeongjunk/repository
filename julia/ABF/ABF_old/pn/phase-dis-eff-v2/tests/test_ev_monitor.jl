using Plots
using Random
using LinearAlgebra
include("../library/abf2_disorder.jl")
include("../library/abf2_ham.jl")
include("../library/abf2_pnscan.jl")

#GIF Save directory
savedir = "/Users/pcs/data/ABF/analysis/nu2-eigvect-sample-gifs/"

function draw_ham(H)
    p = heatmap(real.(Array(H)), c=:thermal)
end

function draw_ham2(H)
    heatmap(abs.(Array(H)), c=:thermal)
end

# function eigenmonitor(H)
function absites(eigvect::AbstractVector)
    n = length(eigvect)
    i_a = collect(2:2:n)
    i_b = collect(1:2:n-1)

    return eigvect[i_b], eigvect[i_a]
end

function ev_monitor(H, str::String; split = true)
    eig = eigen(Hermitian(Array(H)))
    anim = @animate for i in 1:size(H,1)


        pn = round(compute_pn(eig.vectors[:,i]), digits = 2)
        e = round(eig.values[i], digits = 5)

        p1 = scatter(eig.values)
        ylabel!("E")
        if split
            a, b = absites(eig.vectors[:,i])
            p1 = scatter(abs.(a))
            scatter!(p1, abs.(b))
        else
            p1 = scatter(abs.(eig.vectors[:,i]))
        end

        xlabel!("l")
        ylabel!("|ψ_l|")
        title!(str*", PN = $(pn) at E = $(e)")
    end
    return anim

end

function isherm_approx(H)
    return maximum(round.(abs.(H - H'), digits = 13)) == 0
end

#-------------------Parameters-------------------#
θ = 0.25; N = 301; W = 0.00000001; V = 1.
param_str = "th$(θ)-N$(N)-W$(W)-delta$(delta)"
param_str = replace(param_str, "."=>"")

#-------------------Initialize-------------------#
r_on = W*(rand(2N) .- 0.5)
r_off = W*delta*(rand(12N) .- 0.5)


UI = LUT(θ,N)
T = uc_redef(N)
UII = LUT(θ,N)
U = UI*T*UII

#Entangle
H_FD = ham_FD(1,N)
H_FE = U*H_FD*U'

D = spdiagm(0=>r_on)
T_dis = off_phase_dis2(H_FE, r_off)
H_SF = ham_sf_ph(H_FE, D, T_dis, U, W, V)
eig = eigen(Hermitian(Array(H_SF)))

plot(real.(eig.vectors[:,151]))
plot!(imag.(eig.vectors[:,151]))
plot(log.(abs.(eig.vectors[:,151])))


plot!(abs.(eig.vectors[:,151]))

#-----------------------Fully entangled--------------------#

# 1. clean case
str = "clean"
println(str*" is hermitian?"*"$(isherm_approx(H_FE))")
println()
anim = ev_monitor(Hermitian(Array(H_FE)), "clean")
gif(anim, savedir*param_str*"_"*str*".gif", fps = 10)

# 2. onsite disorder case
H_FE_dis = H_FE + W*D
str = "onsite"
println(str*" is hermitian?"*"$(isherm_approx(H_FE_dis))")
println()
anim = ev_monitor(Hermitian(Array(H_FE_dis)), str)
gif(anim, savedir*param_str*"_"*str*".gif", fps = 10)

# 3. onsite and phase disorder case
str = "on_ph"
H_FE_dis = H_FE + D + T_dis
println(str*" is hermitian?"*"$(isherm_approx(H_FE_dis))")
println()
anim = ev_monitor(Hermitian(Array(H_FE_dis)), str)
gif(anim, savedir*param_str*"_"*str*".gif", fps = 10)

# 3. phase disorder only case
str = "ph"
H_FE_dis = H_FE + T_dis
println(str*" is hermitian?"*"$(isherm_approx(H_FE_dis))")
println()
anim = ev_monitor(Hermitian(Array(H_FE_dis)), str)
gif(anim, savedir*param_str*"_"*str*".gif", fps = 10)

#-----------------------Fully detangled--------------------#

# 1. clean case
str = "clean_fd"
H_FD = U'*H_FE*U
println(str*" is hermitian?"*"$(isherm_approx(H_FD))")
println()
anim = ev_monitor(Hermitian(Array(H_FD)), "clean")
gif(anim, savedir*param_str*"_"*str*".gif", fps = 10)

# 2. onsite case
str = "onsite_fd"
H_FD = U'*(H_FE+D)*U
println(str*" is hermitian?"*"$(isherm_approx(H_FD))")
println()
anim = ev_monitor(Hermitian(Array(H_FD)), "clean")
gif(anim, savedir*param_str*"_"*str*".gif", fps = 10)

# 3. onsite + phase case
str = "on_ph_fd"
H_FD = U'*(H_FE+D+T_dis)*U
println(str*" is hermitian?"*"$(isherm_approx(H_FD))")
println()
anim = ev_monitor(Hermitian(Array(H_FD)), "clean")
gif(anim, savedir*param_str*"_"*str*".gif", fps = 10)


# 4. phase only case
str = "ph_fd"
H_FD = U'*(H_FE+D)*U
println(str*" is hermitian?"*"$(isherm_approx(H_FD))")
println()
anim = ev_monitor(Hermitian(Array(H_FD)), "clean")
gif(anim, savedir*param_str*"_"*str*".gif", fps = 10)

#-----------------------Projected detangled--------------------#

# 1. clean case
str = "clean_fd_proj"
H_FD = projection(U'*H_FE*U)
println(str*"is hermitian?"*"$(isherm_approx(H_FD))")
println()
anim = ev_monitor(Hermitian(Array(H_FD)), "clean", split = false)
gif(anim, savedir*param_str*"_"*str*".gif", fps = 10)

# 2. onsite case
str = "onsite_fd_proj"
H_FD = projection(U'*(H_FE+D)*U)/W
println()
println(str*" is hermitian?"*"$(isherm_approx(H_FD))")
anim = ev_monitor(Hermitian(Array(H_FD)), "clean", split = false)
gif(anim, savedir*param_str*"_"*str*".gif", fps = 10)

# 3. onsite + phase case
str = "on_ph_fd_proj"
H_FD = projection(U'*(H_FE+D+T_dis)*U)/W
println()
println(str*" is hermitian?"*"$(isherm_approx(H_FD))")
anim = ev_monitor(Hermitian(Array(H_FD)), "clean", split = false)
gif(anim, savedir*param_str*"_"*str*".gif", fps = 10)


# 4. phase only case
str = "ph_fd_proj"
H_FD = projection(U'*(H_FE+D)*U)
println()
println(str*" is hermitian?"*"$(isherm_approx(H_FD))")
anim = ev_monitor(Hermitian(Array(H_FD)), "clean", split = false)
gif(anim, savedir*param_str*"_"*str*".gif", fps = 10)
