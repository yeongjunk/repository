# Bunch decompositions

using  LinearAlgebra, SparseArrays
using LoopVectorization
export Bunch, bunch!, bunchrank2update!, bunch2!, skewger2!


struct Bunch{T, L<:UnitLowerTriangular{<:T}, D<:Tridiagonal{<:T}, P<:AbstractVector{<:Integer}}
    L::L
    D::D
    p::P
    function Bunch{T,L,D,P}(Lm,Dm,pv) where {T,L,D,P}
        new{T,L,D,P}(Lm, Dm, pv)
    end
end

"""
Finds the index of the absolute maximum of the lower triangular part
"""
function lt_findabsmax(A::AbstractMatrix{T}) where {T <: Real}
    P = size(A, 1)
    Q = size(A, 2)
    m = zero(T);
    ind = (0, 0)
    @inbounds for q in 1:Q, p in q+1:P
        x = abs(A[p,q])
        m >= x && continue
        m = x
        ind = (p, q)
    end 
    return m, ind 
end

"""
Rank 2 update of the skew symmetric matrix A -> A + a*(x*y' - y*x')
"""
@inline function skewger2!(A::AbstractMatrix, a::Real, x::AbstractVector, y::AbstractVector)
    n = size(A, 1)
    @inbounds for j in 1:2:2*div(n, 2)
        xj, yj   = x[j], y[j]
        xj1, yj1 = x[j+1], y[j+1]
        @inbounds for k in j+1:n
            xk, yk    = x[k], y[k]
            A[k, j]   += a*(xk*yj-yk*xj)
            A[k, j+1] += a*(xk*yj1-yk*xj1)
        end
     end
end


function Bunch(L::UnitLowerTriangular{<:T}, D::Tridiagonal{<:T}, p::AbstractVector{<:Integer}) where {T<:Real}
    return Bunch{T, typeof(L), typeof(D), typeof(p)}(L, D, p)
end

function bunch!(A::AbstractMatrix{<:Real}; pivot::Bool=true)
    n = size(A, 1)
    m = size(A, 1)÷2
    rowperm = collect(1:n)
    
    @inbounds for i in 1:2:n-1
        if pivot
            T = view(A, i:n, i:i+1)
            maxval, idx = lt_findabsmax(T)
            if !(maxval <= abs(T[2, 1]))
                q = idx[1] + i-1
                p = (idx[2] == 2) ? i : i+1
                A[q, p] = -A[q, p]
                @inbounds for t in 1:p-1
                    A[p, t], A[q, t] = A[q, t], A[p, t]
                end
                @inbounds for t in p+1:q-1
                    A[t, p], A[q, t] = -A[q, t], -A[t, p]
                end
                @inbounds for t in q+1:n
                    A[t, p], A[t, q] = A[t, q], A[t, p]
                end
                rowperm[p], rowperm[q] = rowperm[q], rowperm[p]
            end
        end
        
        a = A[i+1, i]
        ASinv1 = view(A, i+2:n, i)
        ASinv2 = view(A, i+2:n, i+1)
        B = view(A, i+2:n, i+2:n)

        @inbounds for j in 1:n-i-1
            x, y = ASinv1[j], ASinv2[j]
            ASinv1[j], ASinv2[j] = -y/a, x/a
        end
        # mul!(B, ASinv1, ASinv2', a, 1.0)
        skewger2!(B, a, ASinv1, ASinv2) 
    end
    D = Tridiagonal(A)
    D.dl[2:2:end] .= zero(eltype(D))
    D.du .= -D.dl
    D.d .= zero(eltype(D))

    @views A[diagind(A, -1)][1:2:end] .= zero(eltype(D)) 
    L = UnitLowerTriangular(A);

    return Bunch(L, D, rowperm)
end


function bunch2!(A::AbstractMatrix{<:Real}; pivot=true)
    n = size(A, 1)
    m = size(A, 1)÷2
    rowperm = collect(1:n)
    
    @inbounds for i in 1:2:n-1
        if pivot
            T = view(A, i:n, i:n)
            m, idx = lt_findabsmax(T)

            q = idx[1] + i-1
            p = i
            q1 = idx[2] != 1 ? idx[2] + i-1 : idx[1] + i-1 
            p1 = i+1 

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

            A[q1, p1] = -A[q1, p1]
            @inbounds for t in 1:p1-1
                A[p1, t], A[q1, t] = A[q1, t], A[p1, t]
            end
            @inbounds for t in q1+1:n
                A[t, p1], A[t, q1] = A[t, q1], A[t, p1]
            end
            @inbounds for t in p1+1:q1-1
                A[t, p1], A[q1, t] = -A[q1, t], -A[t, p1]
            end

            rowperm[p],  rowperm[q]  = rowperm[q],  rowperm[p]
            rowperm[p1], rowperm[q1] = rowperm[q1], rowperm[p1]
        end
        
        a = A[i+1, i]
        ASinv1 = view(A, i+2:n, i)
        ASinv2 = view(A, i+2:n, i+1)
        B = view(A, i+2:n, i+2:n)

        @inbounds for j in 1:n-i-1
            x, y = ASinv1[j], ASinv2[j]
            ASinv1[j], ASinv2[j] = -y/a, x/a
        end
        skewger2!(B, a, ASinv1, ASinv2) 

    end
    D = Tridiagonal(A)
    D.dl[2:2:end] .= zero(eltype(D))
    D.du .= -D.dl
    D.d .= zero(eltype(D))

    @views A[diagind(A, -1)][1:2:end] .= zero(eltype(D)) 
    L = UnitLowerTriangular(A);
    
    return Bunch(L, D, rowperm)
end

function bunch!(A, l; pivot=true)
    n = size(A, 1)
    m = size(A, 1)÷2
    rowperm = collect(1:n)
    
    @inbounds for i in 1:2:n-1
        nzlen = min(l^2, n-i-1)
        if pivot
            T = view(A, i:i+nzlen, i:i+1)
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

        @inbounds for j in 1:2:nzlen
            x, y   = ASinv1[j],   ASinv2[j]
            x1, y1 = ASinv1[j+1], ASinv2[j+1]
            @inbounds for k in j+1:nzlen
                z, w = ASinv1[k], ASinv2[k]
                B[k, j]   -= (w*x - z*y)/a
                B[k, j+1] -= (w*x1 - z*y1)/a
            end
        end
        @inbounds for j in 1:nzlen
            x, y = ASinv1[j], ASinv2[j]
            ASinv1[j], ASinv2[j] = -y/a, x/a
        end
    end
    D = Tridiagonal(A)
    D.dl[2:2:end] .= zero(eltype(D))
    D.du .= -D.dl
    D.d .= zero(eltype(D))

    A[diagind(A, -1)] -= D.dl
    L = UnitLowerTriangular(A);
    return Bunch(L, D, rowperm)
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

