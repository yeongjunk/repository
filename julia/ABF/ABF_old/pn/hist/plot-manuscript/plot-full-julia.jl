using Gadfly
using ColorSchemes
using JLD
using Interpolations
import Cairo, Fontconfig
include("../analysis/xi_abf2.jl")

#data save directory
dir_full = "/Users/pcs/data/ABF/rawdata/nu2-tm/nu2-tm-semi-detangled/"
#read all the data and the parameter
ξ, params_weak = ξ_read(dir_full*"/config-xi_all.mat")
(E, W, θI, θII) =params_expand(params_weak)



z = [ξ[:,:,i,i] for i in [2 10 18 26]]

A = reshape(collect(Iterators.product(E,W)), :,1)
x = first.(A)
y = last.(A)

# set_default_plot_size(9cm, 8cm)

latex_fonts = Theme(
    panel_fill= "white",
    major_label_font="Times", major_label_font_size=20pt,
    minor_label_font="Times", minor_label_font_size=20pt,
    key_title_font="Times", key_title_font_size=16pt,
    key_label_font="Times", key_label_font_size=14pt,
    panel_stroke = "black",
    panel_line_width = 2pt,
    plot_padding = [1mm, 1mm, 1mm, 1mm],
    minor_label_color = "black",
    major_label_color = "black",
)


latex_fonts_nonekey = Theme(
    panel_fill= "white",
    major_label_font="Times", major_label_font_size=17pt,
    minor_label_font="Times", minor_label_font_size=17pt,
    key_title_font="Times", key_title_font_size=17pt,
    key_label_font="Times", key_label_font_size=10pt,
    key_position = :none,
    panel_stroke = "black",
    panel_line_width = 2pt,
    plot_padding = [1mm, 1mm, 1mm, 1mm],
    minor_label_color = "black",
    major_label_color = "black",
)


Gadfly.push_theme(latex_fonts)

p1 = plot(x = x, y = y, color = z[1], Geom.rectbin,
    Coord.cartesian(xmin=x[1], xmax=x[end], ymin=y[1], ymax=y[end]),
    Scale.ContinuousColorScale(p -> get(ColorSchemes.afmhot, p)),
    Guide.xlabel("E"),
    Guide.colorkey(title="ξ value"),
    Guide.ylabel("W"),
    )
p2 = plot(x = x, y = y, color = z[2], Geom.rectbin,
    Coord.cartesian(xmin=x[1], xmax=x[end], ymin=y[1], ymax=y[end]),
    Scale.ContinuousColorScale(p -> get(ColorSchemes.afmhot, p)),
    Guide.xlabel("E"),
    Guide.colorkey(title="ξ value"),
    Guide.ylabel("W"),
    )
p3 = plot(x = x, y = y, color = z[3], Geom.rectbin,
    Coord.cartesian(xmin=x[1], xmax=x[end], ymin=y[1], ymax=y[end]),
    Scale.ContinuousColorScale(p -> get(ColorSchemes.afmhot, p)),
    Guide.xlabel("E"),
    Guide.colorkey(title="ξ value"),
    Guide.ylabel("W"),
    )
p4 = plot(x = x, y = y, color = z[4], Geom.rectbin,
    Coord.cartesian(xmin=x[1], xmax=x[end], ymin=y[1], ymax=y[end]),
    Scale.ContinuousColorScale(p -> get(ColorSchemes.afmhot, p)),
    Guide.xlabel("E"),
    Guide.colorkey(title="ξ value"),
    Guide.ylabel("W"),
    )

p_stack = gridstack([p1 p2; p3 p4])


cd()

draw(PNG("full.PNG", 24cm,16cm ), p_stack)
Gadfly.pop_theme()
Gadfly.push_theme(latex_fonts_nonekey)

l1 = layer(x = x, y = y, color = z[1], Geom.rectbin)
l2 = layer(x = x, y = y, color = z[2], Geom.rectbin)
l3 = layer(x = x, y = y, color = z[3], Geom.rectbin)
l4 = layer(x = x, y = y, color = z[4], Geom.rectbin)

maxvalpt = findmax.(z)
maxpt = last.(maxvalpt)

max1 = layer(x = [E[maxpt[1][1]]], y = [W[maxpt[1][2]]], Geom.point)
max2 = layer(x = [E[maxpt[2][1]]], y = [W[maxpt[2][2]]], Geom.point)
max3 = layer(x = [E[maxpt[3][1]]], y = [W[maxpt[3][2]]], Geom.point)
max4 = layer(x = [E[maxpt[4][1]]], y = [W[maxpt[4][2]]], Geom.point)
p1 = plot(max1,l1,
    Coord.cartesian(xmin=x[1], xmax=x[end], ymin=y[1], ymax=y[end]),
    Scale.ContinuousColorScale(p -> get(ColorSchemes.afmhot, p)),
    Guide.xlabel("E"),
    Guide.colorkey(title="ξ value"),
    Guide.ylabel("W")
)
p2 = plot(max2,l2,
    Coord.cartesian(xmin=x[1], xmax=x[end], ymin=y[1], ymax=y[end]),
    Scale.ContinuousColorScale(p -> get(ColorSchemes.afmhot, p)),
    Guide.xlabel("E"),
    Guide.colorkey(title="ξ value"),
    Guide.ylabel("W")
)
p3 = plot(max3,l3,
    Coord.cartesian(xmin=x[1], xmax=x[end], ymin=y[1], ymax=y[end]),
    Scale.ContinuousColorScale(p -> get(ColorSchemes.afmhot, p)),
    Guide.xlabel("E"),
    Guide.colorkey(title="ξ value"),
    Guide.ylabel("W")
)
p4 = plot(max4,l4,
    Coord.cartesian(xmin=x[1], xmax=x[end], ymin=y[1], ymax=y[end]),
    Scale.ContinuousColorScale(p -> get(ColorSchemes.afmhot, p)),
    Guide.xlabel("E"),
    Guide.colorkey(title="ξ value"),
    Guide.ylabel("W")
)
draw(PDF("full-1.PDF", 10cm,8cm,dpi = 500 ), p1)
draw(PDF("full-2.PDF", 10cm,8cm,dpi = 500 ), p2)
draw(PDF("full-3.PDF", 10cm,8cm,dpi = 500 ), p3)
draw(PDF("full-4.PDF", 10cm,8cm,dpi = 500 ), p4)

Gadfly.pop_theme()
println("data saved")
