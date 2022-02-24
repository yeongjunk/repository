# Custom modules
using ABF
using Lattice
using PN

function box_inds(ltc, b)
    L = ltc.N
    U = ltc.U
    box_pts = b
    box_num = U*L ÷ b

    inds = Array{Int64}(undef, box_pts, box_num)
    box_count = 1
    for x in 1:b:(L-b+1), u in 1:U
        for m in 0:b-1
            inds[m + 1, box_count] = index(ltc, (x + m, u))
        end
        box_count += 1
    end
    return inds
end

# @doc"""
# ltc: 2D Lattice
# psi: wavefunction vector on the lattice support
# idx: box indices
# q: GIPR exponent
# """
function compute_box_iprs(ltc, psi, idx::Array{Int64}; q = 2)
    p = Array{Float64}(undef, size(idx, 2), size(psi, 2))
    for i in 1:size(psi, 2)
        for j in 1:size(idx, 2)
            @views p[j, i] = sum(abs2,  psi[idx[:, j], i])
        end
    end
    return compute_iprs2(p, q = q)
end


function compute_box_iprs(ltc, psi, b::Int64; q = 2)
    idx = box_inds(ltc, b)
    p = Array{Float64}(undef, size(idx, 2), size(psi, 2))
    for i in 1:size(psi, 2)
        for j in 1:size(idx, 2)
            @views p[j, i] = sum(abs2, psi[idx[:, j], i])
        end
    end
    return compute_iprs2(p, q = q)
end
