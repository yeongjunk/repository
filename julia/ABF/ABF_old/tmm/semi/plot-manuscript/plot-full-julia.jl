using Plots
using ColorSchemes
using MAT
using LaTeXStrings
using Plots.PlotMeasures

function relative(f, r; sp)
    p = plot!()
    lims = f(p[sp])
    return lims[1] + r * (lims[2]-lims[1])
end
relativex(r; sp::Int=1) = relative(Plots.xlims, r; sp=sp)
relativey(r; sp::Int=1) = relative(Plots.ylims, r; sp=sp)


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


include("../analysis/xi_abf2.jl")

#data save directory
dir_full = "/Users/pcs/data/ABF1D/rawdata/nu2-tm/nu2-tm-semi-detangled/"
#read all the data and the parameter
ξ, params_weak = ξ_read(dir_full*"/config-xi_all.mat")
(E, W, θI, θII) = params_expand(params_weak)


p_a = heatmap(W, E, ξ[:,:,1,1], cb = :right)
annotate!(p_a, sp=1,[(relativex(0.03; sp=1), relativey(0.93; sp=1), text(L"(a)",:white, :left, 22, "computer modern"))])
annotate!(p_a, sp=1,[(relativex(0.12; sp=1), relativey(0.93; sp=1), text(L"\theta = %$(round((θI[1]/pi)pi, digits = 5))",:white, :left, 22, "computer modern"))])


p_b = heatmap(W, E, ξ[:,:,end÷3÷2+1,end÷3÷2+1], cb = :right)
annotate!(p_b, sp=1,[(relativex(0.03; sp=1), relativey(0.93; sp=1), text(L"(b)",:white, :left, 22, "computer modern"))])
annotate!(p_b, sp=1,[(relativex(0.12; sp=1), relativey(0.93; sp=1), text(L"\theta = %$(round(θI[end÷3÷2+1]/pi, digits = 3))\pi",:white, :left, 22, "computer modern"))])

p_c = heatmap(W, E, ξ[:,:,end÷3+1,end÷3+2+1], cb = :right)
annotate!(p_c, sp=1,[(relativex(0.03; sp=1), relativey(0.93; sp=1), text(L"(c)",:black, :left, 22, "computer modern"))])
annotate!(p_c, sp=1,[(relativex(0.12; sp=1), relativey(0.93; sp=1), text(L"\theta = %$(round(θI[end÷3+1]/pi, digits = 3))\pi",:black, :left, 22, "computer modern"))])

p_d = heatmap(W, E, ξ[:,:,end÷2+1,end÷2+1], cb = :right)
xlabel!(L"W")
ylabel!(L"E")
plot!(right_margin = 2mm)
annotate!(p_d, sp=1,[(relativex(0.03; sp=1), relativey(0.93; sp=1), text(L"(d)" , :black, :left, 22, "computer modern"))])
annotate!(p_d, sp=1,[(relativex(0.12; sp=1), relativey(0.93; sp=1), text(L"\theta = %$(round(θI[end÷2+1]/pi, digits = 3))\pi",:black, :left, 22, "computer modern"))])


plot(p_a, p_b, p_c, p_d, size = (1500, 800), bottom_margin = 2mm, left_margin = 5mm)

surface(
  W, E, ξ[:,:,1, 1],
  c=:viridis, legend=:none,
  # vvvvvvvvvvvv series[:extra_kwargs] vvvvvvvvvvvvv
  nx=50, ny=50, display_option=Plots.GR.OPTION_SHADED_MESH,
  grid = true,
  camera = (40, 70)
)
