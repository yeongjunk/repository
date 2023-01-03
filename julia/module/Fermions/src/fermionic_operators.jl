using LinearAlgebra,SparseArrays
using Combinatorics

# when dealing with particle number conserved situation, these operators are used.

function c_adj_c!(r, c, v, t, i, j, basis; sorted = false)
    len_basis = length(basis)
    setbit_checker_j = (one(eltype(basis)) << (j-1))
    setbit_checker_i = (one(eltype(basis)) << (i-1))
    
    for k in eachindex(basis)
        @inbounds if check_setbit(basis[k], j, y = setbit_checker_j)
            basis_k = toggle_bit(basis[k], j)
            sign = count_ones_ranged(basis_k, j-1, 1)
            if !check_setbit(basis_k, i, y = setbit_checker_i)
                basis_k = toggle_bit(basis_k, i)
                sign += count_ones_ranged(basis_k, i-1, 1)
                
                if sorted
                    l = searchsortedfirst(basis, basis_k)
                else
                    l = findfirst(x -> x == basis_k, basis)
                end
                push!(r, l)
                push!(c, k)
                push!(v, (-1)^(sign)*t)
            end
        end
    end
end

function c_adj_c(t::F, i, j, basis; sorted = false) where F
    len_basis = length(basis)
    r = Int64[]
    c = Int64[]
    v = F[]
    c_adj_c!(r, c, v, t, i, j, basis; sorted = sorted)
    
    sparse(r, c, v, len_basis, len_basis)
end

# Quartic operator
function c_adj_c_adj_c_c!(r, c, v, t, i, j, k, l, basis; sorted = false)
    setbit_checker_i = (one(eltype(basis)) << (i-1))
    setbit_checker_j = (one(eltype(basis)) << (j-1))
    setbit_checker_k = (one(eltype(basis)) << (k-1))
    setbit_checker_l = (one(eltype(basis)) << (l-1))

    len_basis = length(basis)
    for n in 1:len_basis
        if check_setbit(basis[n], l, y = setbit_checker_l) 
            sign = 0
            basis_n = basis[n]
            basis_n = toggle_bit(basis_n, l)
            sign   += count_ones_ranged(basis_n, l-1, 1)
            if check_setbit(basis_n, k, y = setbit_checker_k)
                basis_n = toggle_bit(basis_n, k)  
                sign   += count_ones_ranged(basis_n, k-1, 1)
                if !check_setbit(basis_n, j, y = setbit_checker_j)
                    basis_n = toggle_bit(basis_n, j) 
                    sign   += count_ones_ranged(basis_n, j-1, 1)
                    if !check_setbit(basis_n, i, y = setbit_checker_i)
                        basis_n = toggle_bit(basis_n, i)
                        sign   += count_ones_ranged(basis_n, i-1, 1)
                        if sorted
                            m = searchsortedfirst(basis, basis_n)
                        else
                            m = findfirst(x -> x == basis_n, basis)
                        end
                        push!(r, m)
                        push!(c, n)
                        push!(v, (-1)^(sign)*t)
                    end
                end
            end
        end
    end
end


function c_adj_c_adj_c_c(t::F, i, j, k, l, basis; sorted = false) where F
    len_basis = length(basis)
    r = Int64[]
    c = Int64[]
    v = F[]
    c_adj_c_adj_c_c!(r, c, v, t, i, j, k, l, basis, sorted = sorted)
    sparse(r, c, v, len_basis, len_basis)
end

function c_adj!(r, c, v, t, i, basis; sorted = false)
    len_basis = length(basis)
    setbit_checker_i = (one(eltype(basis)) << (i-1))

    for k in 1:len_basis
        if !check_setbit(basis[k], i, y = setbit_checker_i)
            z = set_bit(basis[k], i)

            if sorted
                l = searchsortedfirst(basis, z)
            else
                l = findfirst(x -> x == z, basis)
            end
            push!(r, l)
            push!(c, k)
            push!(v, (-1)^(count_ones_ranged(z, i-1, 1))*t)
        end
    end
end

function c_adj(t::F, i, basis; sorted = false) where F
    r = Int64[];
    c = Int64[];
    v = F[];
    len_basis = length(basis)
    c_adj!(r, c, v, t, i, basis, sorted=sorted)
    sparse(r, c, v, len_basis, len_basis)
end
