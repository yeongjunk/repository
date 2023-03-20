module TransferMatrices
    import LinearAlgebra: UniformScaling, inv, eltype, I
    import SparseArrays: spdiagm, spzeros
    export tm_inv!, tm_inv, tm!, tm 
    using BlockArrays 

    """
    T: reusable transfer matrix, E: energy, H0: onsite matrix, Tp1:n+1 hopping matrix, Tm1: n-1 hopping matrix
    """
    function tm!(T::AbstractMatrix{F}, E, H0::AbstractMatrix{F}, Tp1::AbstractMatrix{F}, Tm1::AbstractMatrix{F}) where F
        M = size(H0, 1)

        inv_Tp1 = inv(Tp1)
        T[1:M, 1:M] .= inv_Tp1*(H0 - UniformScaling(E))
        T[1:M, M+1:2M] .= inv_Tp1*Tm1
    end

    """
    T: reusable transfer matrix, E: energy, H0: onsite matrix, Tp1:n+1 hopping matrix, Tm1: n-1 hopping matrix
    """
    function tm_inv!(T::AbstractMatrix{F}, E, H0::AbstractMatrix{F}, inv_Tp1::AbstractMatrix{F}, Tm1::AbstractMatrix{F}) where F
        M = size(H0, 1)

        T[1:M, 1:M] .= inv_Tp1*(H0 - UniformScaling(E))
        T[1:M, M+1:2M] .= inv_Tp1*Tm1
    end

    function tm(E::Number, H0::AbstractMatrix{F}, Tp1::AbstractMatrix{F}, Tm1::AbstractMatrix{F}) where F
        M = size(H0, 1)

        inv_Tp1 = inv(Tp1)

        A = inv_Tp1*(H0 - UniformScaling(E))
        B = inv_Tp1*Tm1
        C = spdiagm(0 => fill(F(-1), M))
        T = [A B; C spzeros(F, M, M)]

        return T
    end

    function tm_inv(E::Number, H0::AbstractMatrix{F}, inv_Tp1::AbstractMatrix{F}, Tm1::AbstractMatrix{F}) where F
        M = size(H0, 1)

        A = inv_Tp1*(H0 - UniformScaling(E))
        B = inv_Tp1*Tm1
        C = spdiagm(0 => fill(F(-1), M))
        T = [A B; C spzeros(F, M, M)]

        return T
    end
end
