using ABFSym
using LinearAlgebra, SparseArrays
using Lattice

function makesym!(ltc, H, V1, V2)
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

        H[index(ltc, (n + 1, 1)), index(ltc, (n, 3))] = t1*V2
        H[index(ltc, (n + 1, 1)), index(ltc, (n, 4))] = t2*V2
        H[index(ltc, (n + 1, 2)), index(ltc, (n, 3))] = t3*V2
        H[index(ltc, (n + 1, 2)), index(ltc, (n, 4))] = t4*V2

        H[index(ltc, (n + 1, 3)), index(ltc, (n, 1))] = - t1*V2
        H[index(ltc, (n + 1, 3)), index(ltc, (n, 2))] = - t2*V2
        H[index(ltc, (n + 1, 4)), index(ltc, (n, 1))] = - t3*V2
        H[index(ltc, (n + 1, 4)), index(ltc, (n, 2))] = - t4*V2

        H[index(ltc, (n, 3)), index(ltc, (n + 1, 1))] = t1*V2
        H[index(ltc, (n, 4)), index(ltc, (n + 1, 1))] = t2*V2
        H[index(ltc, (n, 3)), index(ltc, (n + 1, 2))] = t3*V2
        H[index(ltc, (n, 4)), index(ltc, (n + 1, 2))] = t4*V2

        H[index(ltc, (n, 1)), index(ltc, (n + 1, 3))] = -t1*V2
        H[index(ltc, (n, 2)), index(ltc, (n + 1, 3))] = -t2*V2
        H[index(ltc, (n, 1)), index(ltc, (n + 1, 4))] = -t3*V2
        H[index(ltc, (n, 2)), index(ltc, (n + 1, 4))] = -t4*V2


        end
end
