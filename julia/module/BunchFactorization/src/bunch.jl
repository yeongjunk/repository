# Bunch decompositions

using  LinearAlgebra, SparseArrays
export Bunch, bunch!, bunchrank2update!, bunch2!, skewrank2update!, skewfindabsmax, skewrowcolperm


struct Bunch{T, L<:UnitLowerTriangular{<:T}, D<:Tridiagonal{<:T}, P<:AbstractVector{<:Integer}}
    L::L
    D::D
    p::P
    function Bunch{T,L,D,P}(Lm,Dm,pv) where {T,L,D,P}
        new{T,L,D,P}(Lm, Dm, pv)
    end
end

function Bunch(L::UnitLowerTriangular{<:T}, D::Tridiagonal{<:T}, p::AbstractVector{<:Integer}) where {T<:Real}
    return Bunch{T, typeof(L), typeof(D), typeof(p)}(L, D, p)
end

"""
Return (value, index) of the absolute maximum element with value at index(Tuple{Int, Int}) by only accessing lower triangular part of a skew symmetric matrix.
"""
function skewfindabsmax(A)
    P = size(A, 1)
    Q = size(A, 2)
    m = zero(eltype(A));
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
Return (value, index) of the maximum element with value at index(Tuple{Int, Int}), by only accessing lower triangular part of a skew symmetric matrix
"""
function skewfindmax(A)
    P = size(A, 1)
    Q = size(A, 2)
    m = zero(eltype(A));
    ind = (0, 0)
    @inbounds for q in 1:Q, p in q+1:P
        z = A[p, q]
        x = abs(z)
        m >= x && continue
        m = x
        if z > 0
            ind = (p, q)
        else
            ind = (q, p)
        end
    end 
    return m, ind 
end

"""
Returns (value, index) of the maximum element with value at index(Tuple{Int, Int}), by only accessing lower triangular part of a skew symmetric matrix
"""
function skewfindmin(A)
    m, ind = skewfindmax(A)
    return m, reverse(ind)
end

"""
Update lower triangular part of a skew symmetric matrix A -> A + a*(x*y' - y*x').
"""
function skewrank2update!(A, a, x, y)
    n = size(A, 1)
    m = 2*div(n, 2)
    @inbounds for j in 1:2:m
        xj, xj1   = a*x[j], a*x[j+1]
        yj, yj1   = a*y[j], a*y[j+1]
        @inbounds for k in j+1:n
            xk, yk    = x[k], y[k]
            A[k, j]   += xk*yj-yk*xj
            A[k, j+1] += xk*yj1-yk*xj1
        end
    end
end

"""
Update lower triangular part of a skew symmetric matrix A -> A + a*(x*y' - y*x'), where x and y are nonzero up to l.
"""
function bandedskewrank2update!(A, l, a, x, y)
    n = size(A, 1)
    @inbounds for j in 1:min(n, l)
        xj, yj   = x[j], y[j]
        @inbounds for k in j+1:min(n, j+l)
            xk, yk    = x[k], y[k]
            A[k, j]   += a*(xk*yj-yk*xj)
        end
     end
end

"""
Update lower triangular part of a skew symmetric matrix A -> A + a*(x*y' - y*x'). This is fast when x and y have many zeros, but slower if x and y are very dense.
"""
function skewrank2update2!(A, a, x, y)
    n = size(A, 1)

    x_inds = findall(!iszero, x)  
    y_inds = findall(!iszero, y)  
    j_inds = union(x_inds, y_inds) 
    @inbounds for j in j_inds 
        xj, yj   = a*x[j], a*y[j]
        @inbounds for k in j_inds 
            if k > j
                xk, yk    = x[k], y[k]
                A[k, j]   += xk*yj-yk*xj
            end
        end
     end
end




"""
Permute q-th and p-th rows and columns of a skew symmetric matrix A, accessing only the lower triangular part.
"""
function skewrowcolperm!(A, q, p)
    if p == q
        return nothing
    end
    p, q = minmax(p, q)
    n = size(A, 1)
    if q > size(A, 1) || p > size(A, 1)
        error("p,q > size(A, 1)")
    end

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
end

"""
pivot can be :partial, :partial2 (true), or anyhting else will be considered as no pivoting
:partial: fast partial pivoting.
:partial2: stronger pivoting, but permutes the indices strongly.
"""
function bunch!(A; pivot=:partial)
    n = size(A, 1)
    m = size(A, 1)÷2
    permidx = collect(1:n)
    
    @inbounds for i in 1:2:n-1
        if pivot == :partial || pivot == :partial2
            if pivot == :partial
                T = view(A, i:n, i:i+1)
            elseif pivot == :partial2 || pivot == true
                T = view(A, i:n, i:n)
            end
         
            m, idx = skewfindmin(T)
            if m != T[2, 1]
                q = idx[1] + i-1
                p = i
                q1 = idx[2] != 1 ? idx[2] + i-1 : idx[1] + i-1 
                p1 = i+1

                skewrowcolperm!(A, q, p)
                skewrowcolperm!(A, q1, p1)

                permidx[p],  permidx[q]  = permidx[q],  permidx[p]
                permidx[p1], permidx[q1] = permidx[q1], permidx[p1]
            end
        end

        a = A[i+1, i]
        x = view(A, i+2:n, i)
        y = view(A, i+2:n, i+1)
        B = view(A, i+2:n, i+2:n)

        @inbounds for j in 1:n-i-1
            xj, yj = x[j], y[j]
            x[j], y[j] = -yj/a, xj/a
        end
        skewrank2update!(B, a, x, y) 
    end

    D = Tridiagonal(A)
    D.dl[2:2:end] .= zero(eltype(D))
    D.du .= -D.dl
    D.d .= zero(eltype(D))
    @views A[diagind(A, -1)][1:2:end] .= zero(eltype(D)) 
    L = UnitLowerTriangular(A);

    return Bunch(L, D, permidx)
end


"""
Bunch decomposition accessing only lower l subdiagonal.
If no pivot, resulting L has only upto l+1 subdiagonal. (Extremely efficient, but unstable.
If pivoting is used, the band structure is destroyed.
"""
function bunch!(A, l; pivot=:partial)
    n = size(A, 1)
    m = size(A, 1)÷2
    permidx = collect(1:n)

    @inbounds for i in 1:2:n-1
        if pivot == :partial || pivot ==:partial2
            nzlen = n-i-1
        else
            nzlen = min(l+1, n-i-1) 
        end

        if pivot == :partial
            T = view(A, i:n, i:i+1)
            m, idx = skewfindabsmax(T)
            if m != abs(T[2, 1])
                q = idx[1] + i-1
                p = (idx[2] == 2) ? i : i+1
                skewrowcolperm!(A, q, p)
                permidx[p], permidx[q] = permidx[q], permidx[p]
            end

        elseif pivot == :partial2 || pivot == true
            T = view(A, i:n, i:n)
            m, idx = skewfindabsmax(T)

            if m != abs(T[2, 1])
                q = idx[1] + i-1
                p = i
                q1 = idx[2] != 1 ? idx[2] + i-1 : idx[1] + i-1 
                p1 = i+1 
                skewrowcolperm!(A, q, p)
                skewrowcolperm!(A, q1, p1)
                permidx[p],  permidx[q]  = permidx[q],  permidx[p]
                permidx[p1], permidx[q1] = permidx[q1], permidx[p1]
            end
        end

        a = A[i+1, i]
        x = view(A, i+2:n, i)
        y = view(A, i+2:n, i+1)
        B = view(A, i+2:n, i+2:n)
        @inbounds for j in 1:nzlen
            xj, yj = x[j], y[j]
            x[j], y[j] = -yj/a, xj/a
        end
        if pivot == :partial || pivot == :partial2
            skewrank2update2!(B, a, x, y) 
        else
            bandedskewrank2update!(B, nzlen, a, x, y) 
        end
    end
    D = Tridiagonal(A)
    D.dl[2:2:end] .= zero(eltype(D))
    D.du .= -D.dl
    D.d .= zero(eltype(D))

    A[diagind(A, -1)] -= D.dl
    L = UnitLowerTriangular(A);
    return Bunch(L, D, permidx)
end

"""
Update the bunch decomposition of A -> A + a(w*v' - v*w'). If your bunch decomposition is pivoted, you have to permute the w and v somehow. 
"""
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

