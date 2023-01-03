using Fermions
using Random, LinearAlgebra, SparseArrays
using BitIntegers
function test_quadratic(; verbose = false)
    flag = true
    for L in 4:10
        for R in 1:30
            basis = get_basis(L)
            i = rand(1:L)
            j = rand(1:L)
            h1 = c_adj(1.0, i, basis, sorted = true)
            h2 = c_adj(1.0, j, basis, sorted = true)'
            h = h1*h2;

            g = c_adj_c(1.0, i, j, basis, sorted = true)
            if verbose 
                println("L = $(L), i = $(i), j = $(j) : ", all(h .== g))
            end
            flag = flag && all(h .== g)
        end
    end
    return flag
end

function test_quartic(; verbose = false, rng = Random.GLOBAL_RNG)
    flag = true
    for L in 5:8
        for r in 1:30
            basis = get_basis(L)
            i, j, k, l = rand(rng, 1:L), rand(rng, 1:L), rand(rng, 1:L), rand(rng, 1:L)

            h1 = c_adj_c(1.0, i, k, basis, sorted = true)
            h2 = c_adj_c(1.0, j, l, basis, sorted = true)

            h = h1*h2;
            g = -c_adj_c_adj_c_c(1.0, i, j, k, l, basis, sorted = true)
            hh = g - h

            if j == k
                g += c_adj_c(1.0, i, l, basis, sorted = true)
            end

            if verbose 
                println("L = $(L), i = $(i), j = $(j), k = $(k), l = $(l) : ", all(h .== g))
            end
            flag = flag && all(h .== g)
        end
    end
    return flag
end

function test_large()
    L = 255
    basis = sort(get_basis(L, 2; T = UInt256));
    r, c, v = Int64[], Int64[], Float64[]
    for i in 1:L
        c_adj_c!(r, c, v, 1.0, i, mod1(i+1, L), basis; sorted = true)
        c_adj_c!(r, c, v, 1.0, mod1(i+1, L), i, basis; sorted = true)
    end
    H = sparse(r, c, v, length(basis), length(basis))
    
    return true
end

function test_quartic_large(; verbose = false, rng = Random.GLOBAL_RNG)
    flag = true
    for L in 60:10:100
        for r in 1:10
            basis = sort(get_basis(L, 2))
            i, j, k, l = rand(rng, 1:L), rand(rng, 1:L), rand(rng, 1:L), rand(rng, 1:L)

            h1 = c_adj_c(1.0, i, k, basis, sorted = true)
            h2 = c_adj_c(1.0, j, l, basis, sorted = true)

            h = h1*h2;
            g = -c_adj_c_adj_c_c(1.0, i, j, k, l, basis, sorted = true)
            hh = g - h

            if j == k
                g += c_adj_c(1.0, i, l, basis, sorted = true)
            end

            if verbose 
                println("L = $(L), i = $(i), j = $(j), k = $(k), l = $(l) : ", all(h .== g))
            end
            flag = flag && all(h .== g)
        end
    end
    return flag
end

function test_bit_operators()
    x_bit = falses(255)
    x_bit[end] = 1
    x = bitvec_to_int(x_bit; T = UInt256)
    
    a = x .== big(2)^big(254)
    b = check_setbit(x, 255)
    c = toggle_bit(x, 255) == 0
    d = check_setbit(x, 255)
    e = count_ones_ranged(x, 255, 1) == 1
    f = get_basis(300, 1; T = UInt512)[end] .== big(2)^299
    return a && b && c && d && e && f
end

@time begin
    a1 = test_quartic()
    a2 = test_quadratic()
    a3 = test_large()
    a4 = test_quartic_large()
    a5 = test_bit_operators()
    if a1 && a2 && a3 && a4 && a5
        println("Test successful")
    else
        println("Test failed: $([a1 a2 a3 a4 a5] .==true))")
    end
end
