module PN

using LinearAlgebra


"""
Return generalized IPR of eigvects with exponent q. If density = true, it calculates IPR from density
"""
function compute_iprs(eigvects; dims = 1, q = 2, density = false)
    if density
        return vec(sum(x -> x^q, eigvects, dims = dims))
    else
        return vec(sum(x -> abs2(x)^q, eigvects, dims = dims))
    end
end
  
"""
Return IPRs of the density p
"""
function compute_pns(eigvects; dims = 1, q = 2, density = false)
    return 1 ./ compute_iprs(eigvects, dims = dims, q = q, density = density)
end

export compute_iprs, compute_pns, eig_pn!
end
