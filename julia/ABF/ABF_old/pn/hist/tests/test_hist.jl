
using Plots
using Random
using StatsBase
A = collect(0.1:0.1:10)
rng = MersenneTwister(1234)
sample(rng, A, 10)
bins = range(1,100,step = 1)

h1 = fit(Histogram, A, bins, closed = :left)

B = collect(50.1:0.1:60)
bins2 = range(2,100,step = 2)
h2 = fit(Histogram, B, bins, closed =:left)
merge!(h1, h2)
plot(h1.weights)


# 2 dimensional
