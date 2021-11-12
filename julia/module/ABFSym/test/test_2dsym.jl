using ABFSym
using LinearAlgebra, SparseArrays
using Lattice
using Plots
using PN
using Statistics



function symplectic_coupling!(ltc, t, I, J, val, site1, site2, u, v, V)
    push!(I, index(ltc, (site1[1], site1[2], u)))
    push!(J, index(ltc, (site2[1], site2[2], v)))
    push!(val, t*V[1, 1])

    push!(I, index(ltc, (site1[1], site1[2], u)))
    push!(J, index(ltc, (site2[1], site2[2], v + 2)))
    push!(val, t*V[1, 2])

    push!(I, index(ltc, (site1[1], site1[2], u + 2)))
    push!(J, index(ltc, (site2[1], site2[2], v)))
    push!(val, t*V[2, 1])

    push!(I, index(ltc, (site1[1], site1[2], u + 2)))
    push!(J, index(ltc, (site2[1], site2[2], v + 2)))
    push!(val,  t*V[2 ,2])
end

#
# function symplectic_coupling!(ltc, H, site1, site2, u,v, V)
#     t = H[index(ltc, (site1[1], site1[2], u)), index(ltc, (site2[1], site2[2], v))]
#     H[index(ltc, (site1[1], site1[2], u)), index(ltc, (site2[1], site2[2], v))] = t*V[1, 1]
#     H[index(ltc, (site1[1], site1[2], u)), index(ltc, (site2[1], site2[2], v + 2))] = t*V[1, 2]
#     H[index(ltc, (site1[1], site1[2], u + 2)), index(ltc, (site2[1], site2[2], v))] = t*V[2, 1]
#     H[index(ltc, (site1[1], site1[2], u + 2)), index(ltc, (site2[1], site2[2], v + 2))] = t*V[2 ,2]
# end

function symp_v(V1, V2, ϕ)
    V = [V1 exp(-im*ϕ)*V2; -exp(im*ϕ)*V2 V1]
    return V
end

function symp_dis(V1, V2, ϕ)
    V = [V1 exp(-im*ϕ)*V2; -exp(im*ϕ)*V2 V1]
    return V
end


function makesym2d(ltc, H, V1, V2)
    t0 = Float64[]
    t1 = Float64[]
    t2 = Float64[]
    t3 = Float64[]
    t4 = Float64[]

    n, m = 1, 1
    for i in 1:2, j in 1:2
        push!(t0, H[index(ltc, (m, n, i)), index(ltc, (m, n, j))] )
        push!(t1, H[index(ltc, (m, n, i)), index(ltc, (m, n + 1, j))])
        push!(t2, H[index(ltc, (m, n, i)), index(ltc, (m + 1, n, j))])
        push!(t3, H[index(ltc, (m, n, i)), index(ltc, (m + 1, n + 1, j))])
        push!(t4, H[index(ltc, (m + 1, n, i)), index(ltc, (m, n + 1, j))])
    end
    I = Int64[];
    J = Int64[];
    val = ComplexF64[];
    Id = LinearAlgebra.I(2)
    for m in 1:ltc.M, n in 1:ltc.N
        for i in 1:2, j in 1:2
            Vx = symp_v(V1, V2, 0)
            Vy = symp_v(V1, V2, pi/2)
            Vd1 = symp_v(V1, V2, 0)
            Vd2 = symp_v(V1, V2, 0)
            symplectic_coupling!(ltc, t0[2(i-1) + (j-1) + 1], I, J, val, (m, n), (m, n), i, j, Id)
            symplectic_coupling!(ltc, t1[2(i-1) + (j-1) + 1], I, J, val, (m, n), (m, n + 1), i, j, Vx)
            symplectic_coupling!(ltc, t1[2(i-1) + (j-1) + 1], I, J, val, (m, n + 1), (m, n), j, i, Vx')
            symplectic_coupling!(ltc, t2[2(i-1) + (j-1) + 1], I, J, val, (m, n), (m + 1, n), i, j, Vy)
            symplectic_coupling!(ltc, t2[2(i-1) + (j-1) + 1], I, J, val, (m + 1, n), (m, n), j, i, Vy')
            symplectic_coupling!(ltc, t3[2(i-1) + (j-1) + 1], I, J, val, (m, n), (m + 1, n + 1), i, j, Vd1)
            symplectic_coupling!(ltc, t3[2(i-1) + (j-1) + 1], I, J, val, (m + 1, n + 1), (m, n), j, i, Vd1')
            symplectic_coupling!(ltc, t4[2(i-1) + (j-1) + 1], I, J, val, (m + 1, n), (m, n + 1), i, j, Vd2)
            symplectic_coupling!(ltc, t4[2(i-1) + (j-1) + 1], I, J, val, (m, n + 1), (m + 1, n), j, i, Vd2')
        end
    end
    return sparse(I, J, val, size(H, 1), size(H, 2))
end

ltc = Lattice2D(10, 10, 4)
H, U = ham_fe(ltc, -1., 1., 0.125)
H = convert.(ComplexF64, H)
H2 = makesym2d(ltc, H, 1.,  1.)
ishermitian(H2)

# display(Matrix(H[1:4, 1:4]))
# display(Matrix(H2[1:4, 1:4]))
#
# display(Matrix(H[1:4, end-3:end]))
# display(Matrix(H2[1:4, end-3:end]))

D = 0.1*spdiagm(0 => rand(size(H, 1)) .- 0.5)
vals, vecs = eigen(Hermitian((Matrix(H2))))
H_fd = round.(Matrix(U'*H2*U), digits = 12)
scatter(vals)
vals, vecs = eigen(Hermitian((Matrix(project(U*(H2 + 0.1D)*U')))))
pn = compute_pns(vecs)
Ec = mean(vals)
mean(pn[findall(x -> 0.9Ec  < x < 1.1Ec, vals)])

scatter(abs.(vecs[:, 2]), ms = 3)

ths = range(0.01, 0.25, length = 5)
