## An attempt to create a clean ABF where Ando's construction works  in 2D

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



function symp_v(V1, V2, ϕ)
    V = [V1 exp(-im*ϕ)*V2; -exp(im*ϕ)*V2 V1]
    return V
end


function makesym2d(ltc, H, V1, V2, ϕ1, ϕ2, ϕ3, ϕ4)
    t0 = ComplexF64[]
    t1 = ComplexF64[]
    t2 = ComplexF64[]
    t3 = ComplexF64[]
    t4 = ComplexF64[]

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
            Vx = symp_v(V1, V2, ϕ1)
            Vy = symp_v(V1, V2, ϕ2)
            Vd1 = symp_v(V1, V2, ϕ3)
            Vd2 = symp_v(V1, V2, ϕ4)
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


function dis(N, W, rng)
    arr = Vector{Float64}(undef, 4N)
    d_arr = W*(rand(rng, 2N) .- 0.5)
    for i in 1:N
        arr[4(i-1) + 1:4(i-1)+1 + 1] = d_arr[2(i-1) + 1:2(i-1) + 1 + 1]
        arr[4(i-1) + 1 + 2:4(i-1) + 1 + 3] = d_arr[2(i-1) + 1:2(i-1) + 1 + 1]
    end
    return arr
end

function ham_fe_sym(ltc, Ea, Eb, θ, V1, V2)
    H_fd = ham_fd(ltc, -2., 0.)
    U1 = LUT(ltc, 0.20)
    R1 = redef1(ltc)
    R2 = redef2(ltc)
    U_first = U1*R1*U1
    U_rest = U1*R2
    H_sd = U_first*H_fd*U_first'
    H_sd_sc = makesym2d(ltc, H_sd, V1, V2)
    H_fe_sc = U_rest*H_sd_sc*U_rest'
    return H_fe_sc, U_rest*U_first
end

println()
    ϕ1 = 0;
    ϕ2 = pi/2;
    ϕ3 = pi/2;
    ϕ4 = pi/2;

    ltc = Lattice2D(10, 10, 4)
    H_fe, U = ham_fe(ltc, -2., 1., 0.125, pi/2, 0)
    H_fe_re = real.(H_fe)
    H_fe_im = imag.(H_fe)

    droptol!(H_fe_re, 1e-12)
    droptol!(H_fe_im, 1e-12)
    H_fe = H_fe_re + H_fe_im*im
    H_fe_sc = makesym2d(ltc, H_fe, 1., 1., ϕ1, ϕ2, ϕ3, ϕ4)
    vals, vecs = eigen(Hermitian((Matrix(H_fe_sc))))
    scatter(vals)

D = 4Diagonal(dis(100, 1., rng))
vals, vecs = eigen(Hermitian(Matrix(H_fe_sc + D)))
scatter(vals)

pn = compute_pns(vecs)
Ec = mean(vals)
mean(pn[findall(x -> 0.9Ec  < x < 1.1Ec, vals)])

scatter(abs.(vecs[:, 2]), ms = 3)

ths = range(0.01, 0.25, length = 5)
