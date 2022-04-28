module TransferMatrices
    import LinearAlgebra: UniformScaling, inv, eltype, I
    export tm!, tm 
    
    """
    T: reusable transfer matrix, E: energy, H0: onsite matrix, Tp1:n+1 hopping matrix, Tm1: n-1 hopping matrix
    """
    function tm!(T::AbstractMatrix, E::Number, H0::AbstractMatrix, Tp1::AbstractMatrix, Tm1::AbstractMatrix, n::Integer)
        inv_Tp1 = inv(Tp1)
        T[1:n, 1:n] = inv_Tp1*(H0 - UniformScaling(E))
        T[1:n, n+1:2n] = inv_Tp1*Tm1
    end

    """
    Create transfer matrix. See tm!
    """
    function tm(E, H0::AbstractMatrix, Tp1::AbstractMatrix, Tm1::AbstractMatrix, n::Integer)
        T = zeros(eltype(H0), 2size(H0, 1), 2size(H0, 1))
        T[n+1:2n, 1:n] .= I(n)
        tm!(T, E, H0, Tp1, Tm1, n)
        return T
    end
end
