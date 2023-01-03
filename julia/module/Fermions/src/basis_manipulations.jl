using Combinatorics

function check_setbit(x, n; y = nothing)
    if y == nothing
        x & (one(typeof(x)) << (n-1)) != zero(typeof(x))
    else
        x & y != zero(typeof(x))
    end
end

function set_bit(x, n)
    x | one(typeof(x)) << (n-1)
end

function toggle_bit(x, n)
    x ⊻ one(typeof(x)) << (n-1)    
end

function count_ones_ranged(x, l, r)
    one_x = one(typeof(x))
    if l >= r
        s = count_ones(x & (((one_x << l) - one_x)⊻(one_x << (r-1) - one_x)))
    else
        s = 0
    end
    return s
end

function bitvec_to_int(x; T = UInt64)
    s = zero(T)
    for i in 1:length(x.chunks)
        s += T(x.chunks[i]) << (64*(i-1))
    end
    return s
end

function get_basis(L; T = Nothing)
    if T == Nothing
        if L <= 64
            T = UInt64
        elseif 64 < L <= 128
            T = UInt128
        else
            error("L > 128 not automatically supported, if you want, set keyword argument T = UInt512 ")
        end
    end
    
    L = T(L)
        
    return 0:((one(T) << L) - one(T)) 
end

function get_basis(L, N; T = Nothing)
    len_basis = binomial(L, N)
    if len_basis > 10^7
        error("The basis size is too big. This is dangerous. Don't do this.")
    end
    
    if T == Nothing
        if L <= 64
            T = UInt64
        elseif 64 < L <= 128
            T = UInt128
        else
            error("L > 128 not automatically supported, if you want, set keyword argument T = UInt512 ")
        end
    end
    
    a = falses(L)
    for i in 1:N
        a[i]=true
    end
    x = multiset_permutations(a, L)

    out = Vector{T}(undef, len_basis)
    for (i, xi) in enumerate(x)
        out[i] = bitvec_to_int(BitVector(xi), T=T)
    end
    
    return out
end
