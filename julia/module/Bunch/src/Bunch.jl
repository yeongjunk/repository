# Bunch decompositions

module Bunch
using  LinearAlgebra
export bunch!, bunchrank2update!, solve_D

function solve_D(D, b)
    x = similar(b)
    for i in 1:2:size(D, 1)
        x[i] = -b[i+1]/D.du[i]
        x[i+1] = b[i]/D.du[i]
    end
    return x
end

function bunch!(A; pivot=true)
    n = size(A, 1)
    m = size(A, 1)÷2
    rowperm = collect(1:n)
    
    @inbounds for i in 1:2:n-3
        if pivot
            T = view(A, i:n, i:i+1)
            _, idx = findmax(abs, T)
            q = idx[1] + i-1
            p = (idx[2] == 2) ? i : i+1
            
            A[q, p] = -A[q, p]
            @inbounds for t in 1:p-1
                A[p, t], A[q, t] = A[q, t], A[p, t]
            end
            @inbounds for t in q+1:n
                A[t, p], A[t, q] = A[t, q], A[t, p]
            end
            @inbounds for t in p+1:q-1
                A[t, p], A[q, t] = -A[q, t], -A[t, p]
            end
            rowperm[p], rowperm[q] = rowperm[q], rowperm[p]
        end
        
        a = A[i+1, i]
        ASinv1 = view(A, i+2:n, i)
        ASinv2 = view(A, i+2:n, i+1)
        B = view(A, i+2:n, i+2:n)

        @inbounds for j in 1:n-i-1
            ASinv1[j], ASinv2[j] = -ASinv2[j]/a, ASinv1[j]/a
            @inbounds @simd for k in 1:j-1
                B[j, k] += a*(ASinv1[j]*ASinv2[k] - ASinv2[j]*ASinv1[k])
            end
        end
    end
    D = Tridiagonal(A)
    D.dl[2:2:end] .= zero(eltype(D))
    D.du .= -D.dl
    D.d .= zero(eltype(D))

    A[diagind(A, -1)] -= D.dl
    L = UnitLowerTriangular(A);
    
    return rowperm, L, D
end

function bunchrank2update!(L, D, a, w, v)
    n = size(L, 1)
    n == 0 && return nothing
    w1, w2, v1, v2 = w[1], w[2], v[1], v[2]
    t = a*(w1*v2 - v1*w2)
    w3   = view(w, 3:n)
    v3   = view(v, 3:n)
    S1   = view(D, 1:2, 1:2)
    l_11 = view(L, 3:n, 1)
    l_12 = view(L, 3:n, 2)
    s1 = S1[2,1]  
    S1[1,2] -= t
    S1[2,1] += t
    s1u = S1[2,1]
    
    x = v1*l_11+v2*l_12-v3
    y = -w1*l_11-w2*l_12+w3
    b = a*s1/s1u
    
    l_11 .= (s1*l_11 + a*(v2*w3 - w2*v3))/s1u
    l_12 .= (s1*l_12 - a*(v1*w3 - w1*v3))/s1u
    
    bunchrank2update!(view(L, 3:n, 3:n), view(D, 3:n, 3:n), b, x, y, D_only = D_only)
end

end

