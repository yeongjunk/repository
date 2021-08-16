using Plots
using LinearAlgebra
using StatsBase
include("../library/abf2_cores.jl")
include("../library/abf2_ham.jl")
include("../library/abf2_pnscan.jl")

savedir = "/Users/pcs/data/ABF/analysis/nu2-pn-phase-dis-new-results/"
BLAS.set_num_threads(1) #Turn off BLAS multi-threading

@doc """
Compute sum g1 over eigenvalue
"""
function shifter(vect)
    val, i_max = findmax(vect)
    return circshift(vect, -i_max+1)
end

function average_ev(θ, N, R)
    H_fd = ham_FD(1,N)
    U = LUT(θ, N)
    T = uc_redef(N)
    U = U*T*U

    H_fe = U*H_fd*U'
    vectsum = zeros(Float64, N)
    for realizations in 1:R
        r_on = rand(2N) .- 0.5
        r_off = rand(12N) .- 0.5
        D_dis = spdiagm(0 => r_on)
        T_dis = off_phase_dis2(H_fe, r_off)

        H_dis = H_fd + T_dis
        H_out = projection(H_dis)

        eig = eigen(Hermitian(Array(H_out)))
        vectsum .+= shifter(abs.(eig.vectors[:, N÷2+1]))
    end

    return vectsum
end

function foo(binnum)
    N = 51
    n = 10000

    #Hamiltonians and unitaries
    Uni = LUT(0.01, N)
    T = uc_redef(N)
    U_full = Uni*T*Uni
    H_fd = ham_FD(1, N)

    #discretized bins
    E_edges = range(-0.5, 0.5, length = binnum + 1)
    PN_mean = zeros(Float64, binnum)
    num = zeros(Float64, binnum)

    pnzero = 0.; # storage for PN(E=0)

    H_fe = (U_full)*(H_fd)*(U_full)'

    for i in 1:n

        r_on = rand(2*N) .- 0.5
        r_off =  rand(12*N) .- 0.5
        D_dis = spdiagm(0 => r_on)
        T_dis = off_phase_dis2(H_fe, r_off)

        W = 1. ;V = 1.

        E, PN =  pn_sf_ph(H_fe, D_dis, T_dis, U_full, W, V);
        E = E/(maximum(E)-minimum(E))
        counts, binned_PN = binned_sum(E, PN, E_edges, count = true)
        num .+= counts
        PN_mean .+= binned_PN
        pnzero += PN[N÷2+1]
    end
    e = midpoints(E_edges)
    pnmean = PN_mean./num
    pnzero = pnzero/n

    return  e, pnmean, pnzero
end

function eliminate_nnn(H)
    for i in 1:size(H,1)
        H[i,mod1(i-2, size(H,1))] = 0.
        H[i,mod1(i+2, size(H,1))] = 0.
    end
    return H
end

function eliminate_nn(H)
    for i in 1:size(H,1)
        H[i,mod1(i-1, size(H,1))] = 0.
        H[i,mod1(i+1, size(H,1))] = 0.
    end
    return H
end

binnum = [21; 201;]
result = foo.(binnum)

p = scatter(result[2][1], result[2][2],label="#E_bins = 21",framestyle = :box)
# plot!(p, result[3][1], result[3][2],label="#E_bins = 41")
# plot!(p, result[4][1], result[4][2],label="#E_bins = 51")

scatter!(p, [0],[result[1][3]], label = "E = 0", annotations = (-0.3, 3, Plots.text("Θ=π/4, Phase disorder only", :left)))
fn = "/Users/pcs/data/ABF/analysis/pd_valley.png"
xlabel!("Midpoints of E_bins")
ylabel!("<PN>")
savefig(p, fn)


pnzero = 0.; # storage for PN(E=0)



    ## Hamiltonians and unitaries
    N = 601;    θ = 10^-3;
    V = 10^1;    W = 10^-4

    Uni = LUT(θ, N)
    T = uc_redef(N)
    U_full = Uni*T*Uni
    H_fd = ham_FD(1, N)

    H_fe = (U_full)*(H_fd)*(U_full)'

    r_on = rand(2*N) .- 0.5
    r_off =  rand(12*N) .- 0.5
    D_dis = spdiagm(0 => r_on)
    T_dis = off_phase_dis2(H_fe, r_off)

    H_dis = H_fe + V*T_dis
    H_dis = U_full'*H_dis*U_full
    H_out = projection(H_dis)
    # H_out = H_dis
    # H_out = eliminate_nnn(H_out)
    eig = eigen(Hermitian(Array(H_out)))
    eig2 = eigen(Hermitian(Array(H_dis)))
    # _,i = findmin(abs.(eig2.values))
##_______
println("Computed")

    # -------------- EVEN AND ODD SITES -------------- #
    i = N÷2 + 1
    p1 = plot(1:2:length(eig.vectors[:,i]), abs.(eig.vectors[1:2:end,i]), label = "odd sites", markersize = 3)
    plot!(p1, 2:2:length(eig.vectors[:,i]), abs.(eig.vectors[2:2:end,i]), label = "even sites", markersize = 3)
    xlabel!("Sites")
    ylabel!("|ψ_j|")
    annotate!(0.5, 0.5, Plots.text("Θ=$(round(θ, digits = 4))π", :left, "CMU Serif"))
    display(p)

    # -------------- EVEN AND ODD SITES (log scale) -------------- #
    p2 = plot(1:2:length(eig.vectors[:,i]), log.(abs.(eig.vectors[1:2:end,i])), label = "odd sites", markersize = 3)
    plot!(p2, 2:2:length(eig.vectors[:,i]), log.(abs.(eig.vectors[2:2:end,i])), label = "even sites", markersize = 3)
    xlabel!("Sites")
    ylabel!("log[|ψ_j]")
    # ------------------------------------------------ #
    p3 = scatter(eig.values, markersize = 1, legend = false)
    xlabel!("#")
    ylabel!("eigenvalues")

    p4 = scatter(eig.values[N÷2-4:N÷2+6], markersize = 2, legend = false)
    xlabel!("#")
    ylabel!("eigenvalues")

    plot(p1, p2, p3, p4)
