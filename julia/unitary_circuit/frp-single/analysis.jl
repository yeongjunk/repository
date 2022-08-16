using Plots


include("lib.jl")

marker = (:circle, 4.5, 1., stroke(-0.5, 1., :black))
line = (:line, :solid, 2.5)
default(
    framestyle = :box,
    size = (600,400),
    # right_margin = [3mm 0mm],
    grid = false,
    minorticks = true,
    legend = (0.1, 0.72),
    fontfamily = "computer modern",
    tickfontsize = 19, 
    guidefontsize = 22, 
    legendfontsize = 18, 
    annotationfontsize = 20, palette = :default)


