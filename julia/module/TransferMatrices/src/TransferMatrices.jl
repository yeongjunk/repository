module TransferMatrices
    import LinearAlgebra: UniformScaling, inv, eltype, I
    export tm!, tm 
    using BlockArrays 

    """
    T: reusable transfer matrix, E: energy, H0: onsite matrix, Tp1:n+1 hopping matrix, Tm1: n-1 hopping matrix
    """
    function tm!(T::AbstractMatrix{F}, E::Number, H0::AbstractMatrix{F}, Tp1::AbstractMatrix{F}, Tm1::AbstractMatrix{F}) where F
        @assert size(T, 2) == size(T, 1) == 2size(H0, 1) == 2size(Tp1, 1) == 2size(Tm1, 1) == 2size(H0, 2) == 2size(Tp1, 2) == 2size(Tm1, 2)
        M = size(H0, 1)

        inv_Tp1 = inv(Tp1)
        T[1:M, 1:M] .= inv_Tp1*(H0 - UniformScaling(E))
        T[1:M, M+1:2M] .= inv_Tp1*Tm1
    end

    function tm(E::Number, H0::AbstractMatrix{F}, Tp1::AbstractMatrix{F}, Tm1::AbstractMatrix{F}) where F
        @assert size(T, 2) == size(T, 1) == 2size(H0, 1) == 2size(Tp1, 1) == 2size(Tm1, 1) == 2size(H0, 2) == 2size(Tp1, 2) == 2size(Tm1, 2)
        M = size(H0, 1)

        inv_Tp1 = inv(Tp1)

        A = inv_Tp1*(H0 - UniformScaling(E))
        B = inv_Tp1*Tm1
        C = spdiagm(0 => fill(F(-1), M))
        T = [A B; C spzeros(F, M, M)]

        return T
    end
end
