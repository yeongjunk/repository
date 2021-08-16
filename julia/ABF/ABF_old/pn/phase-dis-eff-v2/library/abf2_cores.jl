using LinearAlgebra
using SparseArrays

## off-diagonal disorders ##
@doc """
offdiagonal phase disorder ver. 1
"""
function off_phase_dis!(H, r_off)
    I, J, V = findnz(H)
    for i in 1:length(V)
        if I[i] > J[i] #exclude diagonals, lower triangle
            H[I[i],J[i]] = V[i]*exp(im*r_off[i])
            H[J[i],I[i]] = conj(H[I[i],J[i]])
        end
    end
end


@doc """
offdiagonal phase disorder ver. 2
"""
function off_phase_dis2!(H, r_off)
    I, J, V = findnz(H)
    for i in 1:length(V)
        if I[i] > J[i] #exclude diagonals, lower triangle
            H[I[i],J[i]] = V[i]*(1. + im*r_off[i])
            H[J[i],I[i]] = conj(H[I[i],J[i]])
        end
    end
end

@doc """
Unlike off_phase_dis2!, this returns 'phase disorder matrix' only, but
it does not add to the original hamiltonian.
"""
function off_phase_dis2(H, r_off)
    H_copy = copy(H)
    rows = rowvals(H_copy)
    vals = nonzeros(H_copy)

    for j = 1:size(H_copy, 1)
       for i in nzrange(H_copy, j)
          if j > rows[i]
              vals[i] = vals[i]*im*r_off[i]
              H_copy[j,rows[i]] = conj(vals[i])
          elseif j == rows[i]
              vals[i] = 0.
          end
      end
    end
    dropzeros!(H_copy)
    return H_copy
end

## participation numbers ##

@doc """
computate PN of an eigenvector.
"""
function compute_pn(eigvect)
    PN = sum(x -> abs2(x)^2, eigvect)
    return 1/ PN
end

@doc """
computate PN of multiple eigenvectors. Each eigenvectors must be stored in each columns.
Faster than calling compute_pn multiple times.
"""
function compute_pns(eigvects)
    PN = vec(sum(x -> abs2(x)^2, eigvects, dims = 1))
    return 1 ./ PN
end

@doc """
Compute eigenvalues and eigenvectors, and then compute the PN of each eigenvectors.
returns the arrays (eig.values, eig.PN(eig.vectors)) of given matrix H.
"""
function eig_pn(H, vl, vu)
    Hc = copy(H)
    return eig_pn!(Hc, lv, hv)
end

@doc """
Compute eigenvalues and eigenvectors, and then compute the PN of each eigenvectors.
returns the arrays (eig.values, eig.PN(eig.vectors)) of given matrix H.
"""
function eig_pn(H)
    Hc = copy(H)
    return eig_pn!(Hc)
end



@doc """
In-place eig_pn, returns eigenvalues within range (vl, vu)
"""
function eig_pn!(H, vl, vu)
    eig = eigen!(H, vl, vu)
    PN = compute_pns(eig.vectors)
    return eig.values, PN
end

@doc """
Inplace eig_pn
"""
function eig_pn!(H)
    eig = eigen!(H)
    PN = compute_pns(eig.vectors)
    return eig.values, PN
end


#---------------------SCALE FREE, ON & PHASE---------------------#
@doc """
Add onsite & phase disorder -> detangle -> project(scale free, normalize bandwidth)
"""
function ham_sf_ph_on(H_fe, D_dis, T_dis, U, W, V)
    H_dis = H_fe + D_dis + V/W*T_dis
    H_dis .= U'*H_dis*U
    H_out = projection(H_dis)
    return H_out
end

function pn_sf_ph_on(H_fe, D_dis, T_dis, U, W, V)
    H_out = ham_sf_ph_on(H_fe, D_dis, T_dis, U, W, V)
    (E, PN) = eig_pn!(Hermitian(Matrix(H_out)))
    return E, PN
end

function pn_sf_ph_on(H_fe, D_dis, T_dis, U, W, V, vl, vu)
    H_out = ham_sf_ph_on(H_fe, D_dis, T_dis, U, W, V)
    (E, PN) = eig_pn!(Hermitian(Matrix(H_out)),vl, vu)
    return E, PN
end

#---------------------SCALE FREE, ONLY PHASE---------------------#
@doc """
Add phase disorder -> detangle -> project(scale free. Optionally, E can be normalized)
"""
function ham_sf_ph(H_fe, D_dis, T_dis, U, W, V)
    H_dis = H_fe .+ V*T_dis
    H_dis .= U'*H_dis*U
    # H_out = projection(H_dis)
    H_out = H_dis[2:2:end, 2:2:end]
    return H_out
end

function pn_sf_ph(H_fe, D_dis, T_dis, U, W, V)
    H_out = ham_sf_ph(H_fe, D_dis, T_dis, U, W, V)
    (E, PN) = eig_pn!(Hermitian(Matrix(H_out)))

    return E, PN
end

function pn_sf_ph(H_fe, D_dis, T_dis, U, W, V,  vl, vu)
    H_out = ham_sf_ph(H_fe, D_dis, T_dis, U, W, V)
    (E, PN) = eig_pn!(Hermitian(Matrix(H_out)),vl, vu)

    return E, PN
end

#---------------------FULLY DETANGLED---------------------#
@doc """
Add onsite & phsae return detangled  (0.5 < E < 0.5) (E normalized by bandwidth)
"""

function ham_fd_ph_on(H_fe, D_dis, T_dis, U, W, V)
    H_dis = H_fe + W*D_dis + V*T_dis
    H_out = U'*H_dis*U
    return H_out
end

function pn_fd_ph_on(H_fe, D_dis, T_dis, U, W, V)
    H_out = ham_fd_ph_on(H_fe, D_dis, T_dis, U, W, V)
    (E, PN) = eig_pn!(Hermitian(Matrix(H_out)))

    return E, PN
end

function pn_fd_ph_on(H_fe, D_dis, T_dis, U, W, V, vl, vu)
    H_out = ham_fd_ph_on(H_fe, D_dis, T_dis, U, W, V)
    (E, PN) = eig_pn!(Hermitian(Matrix(H_out)), vl, vu)

    return E, PN
end

@doc """
Add onsite & phsae return detangled  (0.5 < E < 0.5) (E normalized by bandwidth)
"""
function ham_fd_ph(H_fe, D_dis, T_dis, U, W, V)
    H_dis = H_fe + V*T_dis
    H_out = U'*H_dis*U
    return H_out
end

function pn_fd_ph(H_fe, D_dis, T_dis, U, W, V)
    H_out = ham_fd_ph(H_fe, D_dis, T_dis, U, W, V)
    (E, PN) = eig_pn(Hermitian(Matrix(H_out)))

    return E, PN
end

function pn_fd_ph(H_fe, D_dis, T_dis, U, W, V, vl, vu)
    H_out = ham_fd_ph(H_fe, D_dis, T_dis, U, W, V)
    (E, PN) = eig_pn(Hermitian(Matrix(H_out)))

    return E, PN
end

#---------------------FULLY ENTANGLED---------------------#
@doc """
onsite & phsae entangled  (0.5 < E < 0.5) (E can be normalized by bandwidth)
"""
function ham_fe_ph_on(H_fe, D_dis, T_dis, U, W, V)
    H_out = H_fe + W*D_dis + V*T_dis
    return H_out
end

function pn_fe_ph_on(H_fe, D_dis, T_dis, U, W, V)
    H_out = ham_fe_ph_on(H_fe, D_dis, T_dis, U, W, V)
    (E, PN) = eig_pn(Hermitian(Matrix(H_out)))

    return E, PN
end

function pn_fe_ph_on(H_fe, D_dis, T_dis, U, W, V, vl, vu)
    H_out = ham_fe_ph_on(H_fe, D_dis, T_dis, U, W, V)
    (E, PN) = eig_pn(Hermitian(Matrix(H_out)))

    return E, PN
end

@doc """
onsite & phsae entangled  (0.5 < E < 0.5) (E can be normalized by bandwidth)
"""
function ham_fe_ph(H_fe, D_dis, T_dis, U, W, V)
    H_out = H_fe + V*T_dis
    return H_out
end

function pn_fe_ph(H_fe, D_dis, T_dis, U, W, V)
    H_out = ham_fe_ph(H_fe, D_dis, T_dis, U, W, V)
    (E, PN) = eig_pn(Hermitian(Matrix(H_out)))

    return E, PN
end

function pn_fe_ph(H_fe, D_dis, T_dis, U, W, V, vl, vu)
    H_out = ham_fe_ph_on(H_fe, D_dis, T_dis, U, W, V)
    (E, PN) = eig_pn(Hermitian(Matrix(H_out)),vl, vu)

    return E, PN
end
