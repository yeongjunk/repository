module Fermions
include("basis_manipulations.jl")
include("fermionic_operators.jl")

export check_setbit, set_bit, toggle_bit, count_ones_ranged, bitvec_to_int, get_basis, c_adj!, c_adj, c_adj_c!, c_adj_c, c_adj_c_adj_c_c!, c_adj_c_adj_c_c

end
