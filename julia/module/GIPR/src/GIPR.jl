# This is used together with the Lattice package. it is box counting method with lattice index defined in Lattice package convention
module GIPR
using LinearAlgebra
using Lattices

export box_indices, box_coarse
"""
Return box indices with box size b. First dimension is the box indices, second dimension is different boxes
"""
function box_indices(ltc::Lattice1D, b::Integer)
    @assert ltc.N%b == 0
    L = ltc.N
    U = ltc.U
    box_pts = U*b 
    box_num = L÷b 
    inds = Array{Int64}(undef, box_pts, box_num)
    box_count=1 
    for x in 1:b:L-b+1
        for n in 0:b-1, u in 1:U 
            inds[U*n+u, box_count] = index(ltc, (x+n, u)) 
        end 
        box_count += 1
    end 
    return inds
end

"""
Return box indices with box size b. First dimension is the box indices, second dimension is different boxes
"""
function box_indices(ltc::Lattice2D, b::Integer)
    @assert ltc.N%b == 0
    L = ltc.N
    U = ltc.U
    box_pts = U*b^2 
    box_num = L^2÷b^2 
    inds = Array{Int64}(undef, box_pts, box_num)
    box_count=1 
    for x in 1:b:L-b+1, y in 1:b:L-b+1
        idx_count = 1
        for m in 0:b-1, n in 0:b-1, u in 1:U 
            inds[idx_count, box_count] = index(ltc, (x+m,y+n, u)) 
            idx_count += 1
        end 
        box_count += 1
    end 
    return inds
end

"""
Return box indices with box size b. First dimension is the box indices, second dimension is different boxes
"""
function box_indices(ltc::Lattice3D, b::Integer)
    @assert ltc.N%b == 0 && ltc.L == ltc.M == ltc.N
    L = ltc.N
    U = ltc.U
    box_pts = U*b^3
    box_num = L^3÷b^3 
    inds = Array{Int64}(undef, box_pts, box_num)
    box_count=1 
    for x in 1:b:L-b+1, y in 1:b:L-b+1, z in 1:b:L-b+1
        idx_count = 1
        for l in 0:b-1, m in 0:b-1, n in 0:b-1, u in 1:U 
            inds[idx_count, box_count] = index(ltc, (x+l,y+m,z+n, u)) 
            idx_count += 1
        end 
        box_count += 1
    end 
    return inds
end

"""
Transform the density vector(or columns vectors) p into coarsed p with box indices.
"""
function box_coarse(p, box_inds)
    @assert size(p, 1) == length(box_inds) 
    p_coarse  = Array{eltype(p)}(undef, size(box_inds, 2), size(p, 2)) 
    for j in 1:size(p, 2), i in 1:size(box_inds, 2)
        p_coarse[i, j] = sum(p[box_inds[:, i], j])
    end
    
    if typeof(p) <: Vector
        return vec(p_coarse)
    else
        return p_coarse
    end
end
end
