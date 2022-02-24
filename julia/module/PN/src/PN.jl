module PN
    using LinearAlgebra

    function compute_iprs(eigvects; dims = 1, q = 2)
        return vec(sum(x -> abs2(x)^q, eigvects, dims = dims))
    end

    function compute_iprs2(eigvects; dims = 1, q = 2)
        return vec(sum(x -> x^q, eigvects, dims = dims))
    end

    function compute_pns(eigvects; dims = 1, q = 2)
        return 1 ./ compute_iprs(eigvects, dims = dims, q = q)
    end

    function compute_pns2(eigvects; dims = 1, q = 2)
        return 1 ./ compute_iprs2(eigvects, dims = dims, q = q)
    end

    function eig_pn!(H)
        E, vecs = eigen!(H)
        PN = compute_pns(vecs)
        return (E, PN)
    end

    export compute_iprs2, compute_iprs, compute_pns, compute_pns2, eig_pn!
end
