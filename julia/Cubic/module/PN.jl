module PN
    using LinearAlgebra
    using DataFrames

    function compute_iprs(eigvects; dims = 1, q = 2.)
        return vec(sum(x -> abs2(x)^q, eigvects, dims = dims))
    end

    function compute_pns(eigvects; dims = 1, q = 2.)
        return 1 ./ compute_iprs(eigvects, dims = dims, q = q)
    end

    function eig_pn!(H; dims = 1, q = 2.)
        E, vecs = eigen!(H)
        PN = compute_pns(vecs,dims = dims, q = q)
        return (E, PN)
    end

    function eig_pn(H; dims = 1, q = 2.)
        E, vecs = eigen(H)
        PN = compute_pns(vecs,dims = dims, q = q)
        return (E, PN)
    end

    export compute_iprs, compute_pns, eig_pn!, eig_pn
end
