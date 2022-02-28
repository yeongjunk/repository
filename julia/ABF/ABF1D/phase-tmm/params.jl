using Parameters

@with_kw struct Params
	E::Float64 = 0.
	θ::Float64 = 0.25
	V::Float64 = 1.
	W::Float64 = 0.
	q::Int64 = 2
	N::Int64 = 3000
	seed::Int64 = 1234
end
