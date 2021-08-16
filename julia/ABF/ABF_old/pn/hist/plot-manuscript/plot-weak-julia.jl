using Gadfly
using ColorSchemes
using MAT

import Cairo, Fontconfig
include("../analysis/xi_abf2.jl")

#data save directory
dir_weak = "/Users/pcs/data/ABF/rawdata/nu2-tm/nu2-tm-semi-detangled-weak-log/"
dir_full = "/Users/pcs/data/ABF/rawdata/nu2-tm/nu2-tm-semi-detangled/"
#read all the data and the parameter
ξ_weak, params_weak = ξ_read(dir_weak*"/config-weak-xi_all.mat")
(E, W, θI, θII) =params_expand(params_weak)

A = reshape(collect(Iterators.product(θI[2:end-1],θII[2:end-1])), :,1)
x = first.(A)
y = last.(A)
z = ξ_weak[1,2,2:end-1,2:end-1]
z = reshape(z, :,1)

# set_default_plot_size(9cm, 8cm)

latex_fonts = Theme(
    panel_fill= "white",
    major_label_font="CMU Serif", major_label_font_size=20pt,
    minor_label_font="CMU Serif", minor_label_font_size=16pt,
    key_title_font="CMU Serif", key_title_font_size=16pt,
    key_label_font="CMU Serif", key_label_font_size=16pt,
    panel_stroke = "black",
    panel_line_width = 2pt,
    minor_label_color = "black",
    major_label_color = "black"
    )

Gadfly.push_theme(latex_fonts)

p = plot(x = x/pi, y = y/pi, color = z, Geom.rectbin,
    Coord.cartesian(xmin=0, xmax=1/2, ymin=0, ymax=1/2),
    Scale.ContinuousColorScale(p -> get(ColorSchemes.afmhot, p)),
    Guide.xlabel("Θ<sub>I</sub>/π"),
    Guide.ylabel("Θ<sub>II</sub>/π"),
    )
cd()

draw(SVG("test2.SVG"), p)

Gadfly.pop_theme()


println("data saved")
