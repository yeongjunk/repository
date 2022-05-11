using Plots
using Random
using LinearAlgebra
using SparseArrays
using Lattice, ABF
using LaTeXStrings

@doc """
Unlike off_phase_dis2!, this returns 'phase disorder matrix' only, but
it does not add to the original hamiltonian.
"""
function off_phase_dis2(H, r_off)
    H_copy = copy(H)
    rows = rowvals(H_copy)
    vals = nonzeros(H_copy)

    for j = 1:size(H_copy, 1)
       for i in nzrange(H_copy, j)
          if j > rows[i]
              vals[i] = vals[i]*im*r_off[i]
              H_copy[j,rows[i]] = conj(vals[i])
          elseif j == rows[i]
              vals[i] = 0.
          end
      end
    end
    dropzeros!(H_copy)
    return H_copy
end


@doc """
transfer matrix of nth supercell...
"""
function tm(E, t_1, t_2, t_m1, t_m2, tt_1, tt_2, tt_m1, tt_m2)
    TM = Array{ComplexF64}(undef, 4, 4)
    H_n_0 = [E -t_m1; -tt_1 E]
    T_n_1 = [t_2 t_1; 0 tt_2]
    T_n_m1 = [t_m2 0; tt_m1 tt_m2]
    TM[1:2, 1:2] = inv(T_n_1)*H_n_0
    TM[1:2, 3:4] = -inv(T_n_1)*T_n_m1
    TM[3:4, 1:2] = I(2)
    TM[3:4, 3:4] .= 0
    return TM
end

function ham_sf_ph(H_fe, D_dis, T_dis, U, W, V)
    H_dis = H_fe + V*T_dis
    H_dis = U'*H_dis*U
    H_out = projection(H_dis)
    return H_out
end

@doc """
Add onsite & phase disorder -> detangle -> project(scale free, normalize bandwidth)
"""
function ham_sf_ph_on(H_fe, D_ph, D_on, U, W, V)
    H_dis = H_fe + V*D_ph + W*D_on
    H_dis = U'*H_dis*U
    H_out = projection(H_dis)
    return H_out
end

function correlation(x)
    cor = Float64[]
    for i in (length(x)÷4+1):(length(x)÷4*3)
        push!(cor, 0.)
    end
    for i in (length(x)÷4+1):(length(x)÷4*3)
        for r in 1:(length(x)÷4)
            cor[r] += abs(x[i+r] - x[i])
        end
    end
    return cor
end
rng = MersenneTwister(5123);

function tm_wave(E, N, θ, δ, q)

# Results
    data1 = Complex{BigFloat}[];
    data2 = Complex{BigFloat}[];
    data3 = Complex{BigFloat}[];
    data4 = Complex{BigFloat}[];

    ltc = Lattice1D(N, 2)
    H_fe, U = ham_fe(ltc, -1., 1., θ);
    H_fe = convert.(ComplexF64, H_fe)

    r_on = rand(rng, 2N) .- 0.5
    r_off = rand(rng, 12N) .- 0.5

    D = spdiagm(0=>r_on)
    T = off_phase_dis2(H_fe, r_off);
    H_sf = Hermitian(project(U'*(T + δ*D)*U));

    ψⱼ = convert.(BigFloat, I(4))
    for j in 2:2:N
        t₁ⱼm1 = H_sf[mod1(j-1, N),mod1(j, N)]
        t₂ⱼm1 = H_sf[mod1(j-1, N),mod1(j+1, N)]
        tm₁ⱼm1 = H_sf[mod1(j-1, N),mod1(j-2, N)]
        tm₂ⱼm1 = H_sf[mod1(j-1, N),mod1(j-3, N)]

        t₁ⱼ = H_sf[mod1(j,N),mod1(j+1, N)]
        t₂ⱼ = H_sf[mod1(j,N),mod1(j+2, N)]
        tm₁ⱼ = H_sf[mod1(j,N),mod1(j-1, N)]
        tm₂ⱼ = H_sf[mod1(j,N),mod1(j-2, N)]

        eⱼ = H_sf[mod1(j,N),mod1(j,N)]

        Tⱼ = tm(E, t₁ⱼ, t₂ⱼ, tm₁ⱼ, tm₂ⱼ, t₁ⱼm1, t₂ⱼm1,  tm₁ⱼm1, tm₂ⱼm1)
        ψⱼ = Tⱼ*ψⱼ

        # QR
        if j%q == 0
            ψⱼ[:,2] = ψⱼ[:,2] .- (ψⱼ[:,1]'*ψⱼ[:,2])/(ψⱼ[:,1]'*ψⱼ[:,1])*ψⱼ[:,1]
            ψⱼ[:,3] = ψⱼ[:,3] .- (ψⱼ[:,1]'*ψⱼ[:,3])/(ψⱼ[:,1]'*ψⱼ[:,1])*ψⱼ[:,1] .- (ψⱼ[:,2]'*ψⱼ[:,3])/(ψⱼ[:,2]'*ψⱼ[:,2])*ψⱼ[:,2]
            ψⱼ[:,4] = ψⱼ[:,4] .- (ψⱼ[:,1]'*ψⱼ[:,4])/(ψⱼ[:,1]'*ψⱼ[:,1])*ψⱼ[:,1] .- (ψⱼ[:,2]'*ψⱼ[:,4])/(ψⱼ[:,2]'*ψⱼ[:,2])*ψⱼ[:,2].- (ψⱼ[:,3]'*ψⱼ[:,4])/(ψⱼ[:,3]'*ψⱼ[:,3])*ψⱼ[:,3]
            push!(data1, ψⱼ[1,1])
            push!(data2, ψⱼ[1,2])
            push!(data3, ψⱼ[1,3])
            push!(data4, ψⱼ[1,4])
        end
    end
    cor1 = correlation(convert.(Float64, log.(abs.(data1))))
    cor2 = correlation(convert.(Float64, log.(abs.(data2))))

    return [data1 data2 data3 data4], [cor1 cor2]
end


R = 20
data, cor2 = tm_wave(0., 3000, 0.25, 0., 4)
cor = cor2
data, _ = tm_wave(0., 3000, 0.125, 0., 2)
p1 = plot(log.(abs.(data[:, 1])), line = line, legend = false)
plot!(p1, log.(abs.(data[:, 2])), line = line)
plot!(p1, log.(abs.(data[:, 3])), line = line)
plot!(p1, log.(abs.(data[:, 4])), line = line)
ylabel!(L"\ln |\psi|")
xlabel!("sites")

p2 = plot(log.(abs.(data[:, 2])), legend = false)
ylabel!(L"\ln |\psi|")
xlabel!("sites")


data, _ = tm_wave(0., 3000, 0.125, 0.001, 2)
p3 = plot(log.(abs.(data[:, 2])), legend = false)
ylabel!(L"\ln |\psi|")
xlabel!("sites")

for r in 1:R
    data, cor2 = tm_wave(0., 3000, 0.25, 0., 2)
    cor += cor2
end
p4 = plot(cor[1:end÷2,2], legend = false, line = line)
ylabel!("Correlation")
xlabel!("sites")

p5 = plot(cor[1:end÷2,2].^2, legend = false, line = line)
ylabel!("Correlation^2")
xlabel!("sites")


# Energy scan
inte = -(1:0.5:5)
E_p = 10 .^inte
E_m = -E_p
E = vec([E_m; reverse(E_p)])
ξ = Float64[]
for i in 1:length(E)
    slope = 0
    for r in 1:10
        data, cor2 = tm_wave(E[i], 5000, 0.25, 0., 4)
        slope += (log(abs(data[end, 2])) - log(abs(data[1, 2])))/5000
    end
    push!(ξ, 1/slope)
end

p6 = plot(E, ξ*10, line = line, label = L"\theta = \pi/4", legend = :topright)
xlabel!(L"E")
ylabel!(L"\xi")
p7 = plot(log.(E[end÷2 + 1:end]), log.(10ξ[end÷2 + 1:end]), line = line, label = L"\theta = \pi/4", legend = :topright)
xlabel!(L"\log_{10} E")
ylabel!(L"\log_{10} \xi")
(log.(ξ[end÷2 + 1]) - log.(ξ[end÷2 + 2]))/(log.(E[end÷2 + 1]) - log(E[end÷2 + 2]))


savedir = "/Users/pcs/data/ABF-sum/1d-sf-pd-tm/"
savefig(p1, savedir*"Lyapunov.pdf")
savefig(p2, savedir*"wave_function_del0.pdf")
savefig(p3, savedir*"wave_function_del0.001.pdf")
savefig(p4, savedir*"cor.pdf")
savefig(p5, savedir*"cor_sq.pdf")
savefig(p6, savedir*"xi_E.pdf")
savefig(p7, savedir*"loglog_xi_E.pdf")
