using LinearAlgebra
using StatsBase
using LaTeXStrings
using Plots
include("library/abf2_cores.jl")
include("library/abf2_ham.jl")
include("library/abf2_pnscan.jl")

savedir = "/Users/pcs/data/ABF/analysis/nu2-pn-phase-dis-new-results/"
BLAS.set_num_threads(1) #Turn off BLAS multi-threading

@doc """
Compute sum g1 over eigenvalue
"""
function shifter(vect)
    val, i_max = findmax(vect)
    return circshift(vect, -i_max+1)
end

function bin_label(x, x_edges)
    j = 1 #array index for x
    x_lbl = zeros(Int64, length(x))

    while x[j] < x_edges[1]
        j += 1
    end
    for i in 2:length(x_edges)
        while (j <= length(x)) && (x_edges[i-1] <= x[j] < x_edges[i])
            x_lbl[j] = i-1
            j += 1
        end
    end
    return x_lbl
end

function average_ev(θ, N, R)
    E_binnum = 61
    E_edges = range(-0.5, 0.5, length = E_binnum+1)
    H_fd = ham_FD(1,N)
    U = LUT(θ, N)
    T = uc_redef(N)
    U = U*T*U

    H_fe = U*H_fd*U'
    z_vectsum = zeros(Float64, N)
    vectsum = zeros(Float64, N, E_binnum)

    for realizations in 1:R
        r_on = rand(2N) .- 0.5
        r_off = rand(12N) .- 0.5

        D_dis = spdiagm(0 => r_on)
        T_dis = off_phase_dis2(H_fe, r_off)
        H_dis = H_fe + T_dis
        H_fd_dis = U'*H_dis*U

        H_out = projection(H_fd_dis)

        eig = eigen(Hermitian(Array(H_out)))
        eig.values .= eig.values

        lbl = bin_label(eig.values, E_edges)

        z_vectsum .+= shifter(log.(abs.(eig.vectors[:, N÷2+5])))
        # for i in 1:E_binnum
        #     idx = findall(x-> x == i, lbl)
        #     for j in 1:length(idx)
        #         vectsum[:,i] .+= shifter(log.(abs.(eig.vectors[:,j])))
        #     end
        # end
    end

    return z_vectsum/R
end

function cor(x)
    y = log.(abs.(x))
    out_r = zeros(Float64, length(y)÷2)
    for r in 1:length(y)÷2
        for i in 1:length(y)÷2
            out_r[r] += abs(y[i] - y[i+r])
        end
    end
    return out_r/length(y)
end



θ = 0.25;N = 401; R = 500
z_vectsum = average_ev(θ, N, R);

x = 1:length(z_vectsum)
p1 =  plot(x, z_vectsum, label = L"g(l)",framestyle = :box,  dpi = 150, color = "red")
xlabel!(L"l")
ylabel!(L"g(l)")
savefig(p1, savedir*"gle.png")

p2 = plot(log10.(x),log10.(abs.(z_vectsum)), legend = false, dpi = 150,framestyle = :box, color = "red", markersize = 3)
xlabel!(L"\log_{10}{l}")
ylabel!(L"\log_{10}{g(l)}")
savefig(p2, savedir*"gl_sqe.png")

p3 = plot(p1, p2, layout = 2, size = (600, 300))
savefig(p3, savedir*"glglsqe.png")

## Test correlation function


function bar()
    N = 201
    x = range(-20,20, length = N)
    y = zeros(Float64, N÷2)
    for realizations in 1:2000
        r = 20*rand()
        x = range(-20+r,20+r, length = N)
        x = exp.(abs.(x))
        x .= normalize(x)
        y .+= cor(x)
    end
    return y/2000
end
