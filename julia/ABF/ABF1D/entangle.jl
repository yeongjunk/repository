using Random
using Plots
using LinearAlgebra, SparseArrays

push!(LOAD_PATH, "/Users/pcs/codes/chain/ABF/module")
using ABF
using Lattice
using PN


## Draw lattice
function visualize(ltc::Lattice1D, H)
    vts = site.(Ref(ltc), 1:size(H,1))
    uc = getindex.(vts, 2)
    p = plot(axis = nothing, palette = :tab10, size = (500, 200))

    for i in 1:length(vts), j in 1:length(vts)
        if round(abs(H[i,j]), digits = 12) > 0 && [vts[i][1] vts[j][1]] != [1 ltc.N] && [vts[i][1] vts[j][1]] != [ltc.N 1]
            plot!(p, [vts[i]; vts[j]], color = "black", label = :none, lw = 2)
        end
    end
    scatter!(p, vts, ms = 8, msw = 2, c = uc, grid = false, legend = false, label = :none)

    return p
end
## Parameters
θ = 0.23
W = 1.
Ea = -2.
Eb = 0.
N = 10
## Enangle

l = Lattice1D(N, 2)
l_p = Lattice1D(N, 1)

H, U = ham_fe(l, Ea, Eb, θ)
# D = Array{Float64}(undef, 2N)
# D1 =  W*(rand(size(H, 1)÷2, 1) .- 0.5)
# D2 = -D1
# D[1:2:end-1] = D1
# D[2:2:end] = D2

D = rand(size(H, 1)) .- 0.5
H = H + spdiagm(0 => D) # disorder
H = U'*H*U
droptol!(H, 1E-12) # Detangle and project

H = project(H)

# Can we project in entangled hamiltonian?
H, U = ham_fe(l, 0., 1., θ)
H .+= spdiagm(0 => rand(size(H, 1)) .- 0.5)
P_a = ham_fd(l, 0., 1.)
P_a_fe = U*P_a*U'

H_prj = P_a_fe*H*P_a_fe'
droptol!(H, 1E-12)

E, psi = eigen(Hermitian(Matrix(H)))
E_prj, psi_prj = eigen(Hermitian(Matrix(H_prj)))

E_prj = E_prj[findall(a -> abs(a) > 1E-14, E_prj)]
scatter(E, legend = false)
scatter!(E_prj)
