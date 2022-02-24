using Random
using Plots
using LinearAlgebra, SparseArrays

push!(LOAD_PATH, "/Users/pcs/codes/chain/ABF/module")
using ABF2D
using Lattice
using PN

function visualize(ltc::Lattice2D, H)
    vts = site.(Ref(ltc), 1:size(H,1))
    uc = getindex.(vts, 3)
    p = plot(axis = nothing, palette = :tab10)

    for i in 1:length(vts), j in 1:length(vts)
        if round(abs(H[i,j]), digits = 12) > 0 && [vts[i][1] vts[j][1]] != [1 ltc.M] && [vts[i][1] vts[j][1]] != [ltc.M 1]  && ([vts[i][2] vts[j][2]] != [1 ltc.N]) && ([vts[i][2] vts[j][2]] != [ltc.N 1])
            plot!(p, [vts[i]; vts[j]], color = "black", label = :none, lw = 2)
        end
    end
    scatter!(p, vts, ms = 8, msw = 2, c = uc, grid = false, legend = false, label = :none)

    return p
end

θ = 0.125
W = 1.
l = Lattice2D(5, 5, 2)
H_fd = ham_fd(l, -1, 1)
U = LUT(l, θ)
T1 = redef1(l)
T2 = redef2(l)
H_sd1 = U*H_fd*U'
H_sd2 = U*T1*H_sd1*T1'*U'
H_fe = U*T2*H_sd2*T2'*U'

H_dis = H_fe .+ spdiagm(0 => W*(rand(size(H_fe,1)) .- 0.5))
l_p = Lattice2D(5,5,1)

H_p = project(U'*T1'*U'*T2'*U'*H_dis*U*T2*U*T1*U)


vals, vecs = eigen!(Hermitian(Array(H_p)))
scatter(vecs[:,end÷10])

H_dis_sd2 = T2'*U'*H_dis*U*T2

H_dis_sd1 = T1'*U'*T2'*U'*H_dis*U*T2*U*T1
H_dis_fd = U'*H_dis_sd1*U

p1 = visualize(l, H_fd);
plot!(p1, camera = (20, 87), dpi = 150)
xlabel!("x")
ylabel!("y")

p2 = visualize(l, H_sd1);
plot!(p2, camera = (20, 87),dpi = 150)
xlabel!("x")
ylabel!("y")

p3 = visualize(l, H_sd2);
plot!(p3, camera = (20, 87), dpi = 150)
xlabel!("x")
ylabel!("y")

p4 = visualize(l, H_fe);
plot!(p4, camera = (20, 87), dpi = 150)
xlabel!("x")
ylabel!("y")

p5 = visualize(l_p, H_p);
plot!(p5, camera = (20, 87), dpi = 150)
xlabel!("x")
ylabel!("y")

p6 = visualize(l, H_dis_sd2)
plot!(p6, camera = (20, 87), dpi = 150)

p7 = visualize(l, H_dis_sd1)
plot!(p7, camera = (20, 87), dpi = 150)

p8 = visualize(l, H_dis_fd)
plot!(p8, camera = (20, 87), dpi = 150)


savedir = "/Users/pcs/data/ABF2D/entangling/"

savefig(p1, savedir*"fd.png")
savefig(p2, savedir*"sd1.png")
savefig(p3, savedir*"sd2.png")
savefig(p4, savedir*"fe.png")
savefig(p5, savedir*"sf.png")
