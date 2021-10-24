using ABFSym
using LinearAlgebra, SparseArrays
using Lattice

function makesym!(ltc, H, V1, V2, ϕ)
    n = 1
    t1 = H[index(ltc, (n + 1, 1)), index(ltc, (n, 1))]
    t2 = H[index(ltc, (n + 1, 1)), index(ltc, (n, 2))]
    t3 = H[index(ltc, (n + 1, 2)), index(ltc, (n, 1))]
    t4 = H[index(ltc, (n + 1, 2)), index(ltc, (n, 2))]
    for n in 1:ltc.N
        H[index(ltc, (n + 1, 1)), index(ltc, (n, 1))] = t1*V1
        H[index(ltc, (n + 1, 1)), index(ltc, (n, 2))] = t2*V1
        H[index(ltc, (n + 1, 2)), index(ltc, (n, 1))] = t3*V1
        H[index(ltc, (n + 1, 2)), index(ltc, (n, 2))] = t4*V1

        H[index(ltc, (n + 1, 3)), index(ltc, (n, 3))] = t1*V1
        H[index(ltc, (n + 1, 3)), index(ltc, (n, 4))] = t2*V1
        H[index(ltc, (n + 1, 4)), index(ltc, (n, 3))] = t3*V1
        H[index(ltc, (n + 1, 4)), index(ltc, (n, 4))] = t4*V1

        H[index(ltc, (n, 1)), index(ltc, (n + 1, 1))] = t1*V1
        H[index(ltc, (n, 2)), index(ltc, (n + 1, 1))] = t2*V1
        H[index(ltc, (n, 1)), index(ltc, (n + 1, 2))] = t3*V1
        H[index(ltc, (n, 2)), index(ltc, (n + 1, 2))] = t4*V1

        H[index(ltc, (n, 3)), index(ltc, (n + 1, 3))] = t1*V1
        H[index(ltc, (n, 4)), index(ltc, (n + 1, 3))] = t2*V1
        H[index(ltc, (n, 3)), index(ltc, (n + 1, 4))] = t3*V1
        H[index(ltc, (n, 4)), index(ltc, (n + 1, 4))] = t4*V1

        H[index(ltc, (n + 1, 1)), index(ltc, (n, 3))] = t1*exp(-im*ϕ)*V2
        H[index(ltc, (n + 1, 1)), index(ltc, (n, 4))] = t2*exp(-im*ϕ)*V2
        H[index(ltc, (n + 1, 2)), index(ltc, (n, 3))] = t3*exp(-im*ϕ)*V2
        H[index(ltc, (n + 1, 2)), index(ltc, (n, 4))] = t4*exp(-im*ϕ)*V2

        H[index(ltc, (n + 1, 3)), index(ltc, (n, 1))] = - t1*exp(im*ϕ)*V2
        H[index(ltc, (n + 1, 3)), index(ltc, (n, 2))] = - t2*exp(im*ϕ)*V2
        H[index(ltc, (n + 1, 4)), index(ltc, (n, 1))] = - t3*exp(im*ϕ)*V2
        H[index(ltc, (n + 1, 4)), index(ltc, (n, 2))] = - t4*exp(im*ϕ)*V2

        H[index(ltc, (n, 3)), index(ltc, (n + 1, 1))] = t1*exp(im*ϕ)*V2
        H[index(ltc, (n, 4)), index(ltc, (n + 1, 1))] = t2*exp(im*ϕ)*V2
        H[index(ltc, (n, 3)), index(ltc, (n + 1, 2))] = t3*exp(im*ϕ)*V2
        H[index(ltc, (n, 4)), index(ltc, (n + 1, 2))] = t4*exp(im*ϕ)*V2

        H[index(ltc, (n, 1)), index(ltc, (n + 1, 3))] = -t1*exp(-im*ϕ)*V2
        H[index(ltc, (n, 2)), index(ltc, (n + 1, 3))] = -t2*exp(-im*ϕ)*V2
        H[index(ltc, (n, 1)), index(ltc, (n + 1, 4))] = -t3*exp(-im*ϕ)*V2
        H[index(ltc, (n, 2)), index(ltc, (n + 1, 4))] = -t4*exp(-im*ϕ)*V2
        end
end

ltc = Lattice1D(10, 4)
H, U = ham_fe(ltc, -1, 1, 0.25)
H_fd0 = U'*H*U
H = convert.(ComplexF64, H)
makesym!(ltc, H, 1/sqrt(2), 1/sqrt(2), 0)
ishermitian(H)
H_fd = U'*H*U
ishermitian(H_fd)
droptol!(H_fd, 1E-12)
H_fd = Matrix(H_fd)
vals, vecs = eigen(H_fd)
scatter(vals)
U2 = [0. 0.859197; -0.859197 0.]
H2 = H_fd[1:2,1:2]

vals2, vecs2 = eigen(Matrix(Hermitian(H_fd[1:4, 1:4])))
vecs2 = vecs2[:, [1,3,2,4]]
H_diag = vecs2'*H_fd[1:4, 1:4]*vecs2
round.(H_diag, digits = 12)

vecs2'
display(round.(vecs2, digits = 12))
