push!(LOAD_PATH, "$(pwd())/test_module")
using Hello

Hello.hi()

struct OrderedPair
      x::Real
      y::Real
      OrderedPair(x,y) = x > y ? error("out of order") : new(x,y)
end

p = OrderedPair(1,2)
