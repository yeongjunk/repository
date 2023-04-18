export solve_D!, solve_D!, solve_L!, solve_L, solve_U!, solve_U, solve!, solve

function solve_D(D::Tridiagonal, b)
    x = similar(b)
    for i in 1:2:size(D, 1)
        x[i] = -b[i+1]/D.du[i]
        x[i+1] = b[i]/D.du[i]
    end
    return x
end
function solve_D!(x, D, b)
    for i in 1:2:size(D, 1)
        x[i] = -b[i+1]/D.du[i]
        x[i+1] = b[i]/D.du[i]
    end
    x
end

function solve_L!(x, L::UnitLowerTriangular{F}, b::AbstractVector{F}) where F
    n = length(b)
    @inbounds for i in 1:n
        sum = zero(eltype(b))
        @inbounds for j in 1:i-1
            sum += L[i, j]*x[j]
        end
        x[i] = b[i] - sum
    end
    return x
end

function solve_U!(x, U::UnitUpperTriangular{F}, b::AbstractVector{F}) where F
    n = length(b)

    for i in n:-1:1
        sum = zero(eltype(b))
        for j in n:-1:i+1
            sum += U[i, j]*x[j]
        end
        
        x[i] = b[i] - sum
    end
    return x
end



function solve!(x, F::Bunch, b)
    P = sparse(I(length(b)))[:, F.p]
    x .= P'b
    x .= solve_L!(x , F.L, x)
    x .= solve_D(F.D, x)
    x .= solve_U!(x , F.L', x)
    x .= P*x
    return x
end

function solve(F::Bunch, b)
    x = similar(b) 
    solve!(x, F, b)
    return x
end


function solve_L!(x, l::Integer, L::UnitLowerTriangular{F},  b::AbstractVector{F}) where F
    n = length(b)
    @inbounds for i in 1:n
        sum = zero(eltype(b))
        ii = max(i-l, 1)
        @inbounds for j in ii:i-1
            sum += L[i, j]*x[j]
        end
        x[i] = b[i] - sum
    end
    return x
end

function solve_U!(x, l::Integer, U::UnitUpperTriangular{F},  b::AbstractVector{F}) where F
    n = length(b)

    for i in n:-1:1
        sum = zero(eltype(b))
        ii = min(i+l, n)
        for j in ii:-1:i+1
            sum += U[i, j]*x[j]
        end
        
        x[i] = b[i] - sum
    end
    return x
end

function solve!(x, l::Integer, F::Bunch, b)
    P = sparse(I(length(b)))[:, F.p]
    x .= P'b
    x .= solve_L!(x,l,F.L,x)
    x .= solve_D(F.D,x)
    x .= solve_U!(x,l,F.L',x)
    x .= P*x
    return x
end

function solve(F::Bunch, l::Integer, b)
    x = similar(b) 
    solve!(x,l,F,b)
    return x
end

