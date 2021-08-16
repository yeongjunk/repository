using Gadfly
import Cairo, Fontconfig
using CSV
using DataFrames
using ColorSchemes
using LinearAlgebra
topdir = "/Users/pcs/data/ABF/rawdata/nu2-pn/phase-dis-eff2/"
del0dir = "/Users/pcs/data/ABF/rawdata/nu2-pn/phase-dis-eff2/del0/"
W0dir = "/Users/pcs/data/ABF/rawdata/nu2-pn/phase-dis-eff2/W0/"


savedir = "/Users/pcs/data/ABF/analysis/nu2-pn-phase-dis-eff2/"
fn = "pn_dataframe.csv"
df = DataFrame(CSV.File(topdir*fn))
df.E = round.(df.E, digits = 12)
dfe0 = df[df.E .== 0, :] # at E = 0
dfth25 = df[df.Theta .==0.25,:] # at Θ = 0.25π

#-------df δ = 0 -------#
dfdel0 = DataFrame(CSV.File(del0dir*fn))
dfdel0.E = round.(dfdel0.E, digits = 12)

dfdel0e0 = dfdel0[dfdel0.E .== 0, :]

#-------df W = 0 --------#
dfw0 = DataFrame(CSV.File(W0dir*fn))
dfw0.E = round.(dfw0.E, digits = 12)
dfw0e0 = dfw0[dfw0.E .== 0, :]
replace!(dfw0e0.L10_Delta, 2.0 => 5.0)


println("Data loaded")

ft = "CMU Serif";
theme = Theme(
    panel_fill= "white",
    #-----FONT TYPE-----#
    major_label_font= ft,
    minor_label_font= ft,
    key_title_font  = ft,
    key_label_font  = ft,
    #-----FONT SIZE-----#
    major_label_font_size=18pt,
    minor_label_font_size=18pt,
    key_title_font_size=16pt,
    key_label_font_size=16pt,
    #-----FONT COLOR-----#
    minor_label_color = "black",
    major_label_color = "black",
    key_title_color  =  "black",
    key_label_color  =  "black",
    #-----PANELS & GRID----#
    panel_stroke    =   "black",
    panel_line_width = 2pt,
    grid_line_width = 1pt,
    grid_line_style = :solid,
    plot_padding = [1mm, 5mm, 1mm, 1mm],

)
theme_big = Theme(
    panel_fill= "white",
    #-----FONT TYPE-----#
    major_label_font= ft,
    minor_label_font= ft,
    key_title_font  = ft,
    key_label_font  = ft,
    #-----FONT SIZE-----#
    major_label_font_size=11pt,
    minor_label_font_size=11pt,
    key_title_font_size=10pt,
    key_label_font_size=10pt,
    #-----FONT COLOR-----#
    minor_label_color = "black",
    major_label_color = "black",
    key_title_color  =  "black",
    key_label_color  =  "black",
    panel_stroke    =   "black",
    panel_line_width = 1pt ,
    plot_padding = [1mm, 5mm, 1mm, 1mm],
    grid_line_width = 0.5pt,
    grid_line_style = :solid,
)




Gadfly.with_theme(theme) do
#------------------------------ PN vs THETA ------------------------------#
    cscheme = ColorSchemes.jet1
    xticks = collect(0.00:0.05:0.25)
    yticks = 1:2:9
    l1 = layer(dfe0[dfe0.L10_Delta .> -1.2, :], x = :Theta, y = :PN_mean, color=:L10_Delta, Geom.line, order = 1)
    l2 = layer(dfdel0e0, x = :Theta, y = :PN_mean, color=:Delta, Geom.point, order = 2)
    l3 = layer(dfw0e0[dfw0e0.L10_Delta .== 5.0, :], x = :Theta, y = :PN_mean, color=:L10_Delta, Geom.point, order = 3)
    p1 = plot(l1,l2,l3,
        Guide.xticks(ticks = xticks),
        Guide.yticks(ticks = yticks),
        Guide.xlabel("Θ/π"),
        Guide.ylabel("PN"),
        style(line_width=0.6mm),
        Scale.ContinuousColorScale(p -> get(cscheme, p)),
        Guide.colorkey(title="log<sub>10</sub>[δ]")
    );
    display(p1)

#------------------------------ PN vs DELTA ------------------------------#
    xticks = -2:1:5
    yticks = 1:3:10
    p2 = plot(dfe0, x =:L10_Delta, y =:PN_mean, color=:Theta, Geom.line,
        Guide.xticks(ticks = xticks),
        Guide.yticks(ticks = yticks),
        Scale.ContinuousColorScale(p -> get(cscheme, p)),
        Guide.xlabel("log<sub>10</sub>[δ]"),
        Guide.ylabel("PN"),
        Guide.colorkey(title="Θ/π"),
        style(line_width=0.6mm)
    );
    display(p2)

#------------------------------ PN vs THETA, DELTA ------------------------------#
    xticks = collect(0.00:0.05:0.25)
    xticks[1] = 0.01
    yticks = -2:1:5
    p3 = plot(dfe0, color = :PN_mean, y = :L10_Delta, x = :Theta, Geom.rectbin,
        Scale.ContinuousColorScale(p -> get(ColorSchemes.viridis, p)),
        Guide.yticks(ticks=yticks),
        Guide.xticks(ticks=xticks),
        Guide.xlabel("Θ/π"),
        Guide.ylabel("log<sub>10</sub>[δ]"),
        Guide.colorkey(title="PN"),
        Theme(theme_big)
        )
    display(p3)

#------------------------------ PN vs E at Θ=0.25π ------------------------------#
    xticks = -0.5:0.2:0.5
    yticks = 1:3:10
    p4 = plot(dfth25, x =:E, y =:PN_mean, color=:L10_Delta, Geom.line,
        Scale.ContinuousColorScale(p -> get(cscheme, p)),
        Guide.yticks(ticks=yticks),
        Guide.xticks(ticks=xticks),
        Guide.xlabel("E"),
        Guide.ylabel("PN"),
        Guide.colorkey(title="log<sub>10</sub>[δ]"),
        style(line_width=0.6mm)
        )
    display(p4)

#------------------------------ PN vs DELTA at Θ=0.25π------------------------------#
    xticks = -2:1:5
    yticks = 1:3:10
    p5 = plot(dfth25[-0.5 .< dfth25.E .< 0.0, :], x =:L10_Delta, y =:PN_mean, color=:E, Geom.line,
        Scale.ContinuousColorScale(p -> get(cscheme, p)),
        Guide.yticks(ticks=yticks),
        Guide.xticks(ticks=xticks),
        Guide.xlabel("log<sub>10</sub>[δ]"),
        Guide.ylabel("PN"),
        Guide.colorkey(title="E"),
        style(line_width=0.6mm)
        )
    display(p5)
#------------------------------ PN vs E, DELTA at Θ=0.25π ------------------------------#
    xticks = -0.5:0.2:0.5
    yticks = -2:1:5
    p6 = plot(dfth25, x =:E, y =:L10_Delta, color=:PN_mean, Geom.rectbin,
        Scale.ContinuousColorScale(p -> get(cscheme, p)),
        Guide.yticks(ticks=yticks),
        Guide.xticks(ticks=xticks),
        Guide.xlabel("E"),
        Guide.ylabel("log<sub>10</sub>[δ]"),
        Guide.colorkey(title="PN")
        )
    display(p6)

    draw(PDF(savedir*"p1.pdf"), p1)
    draw(PDF(savedir*"p2.pdf"), p2)
    draw(PDF(savedir*"p3.pdf"), p3)
    draw(PDF(savedir*"p4.pdf"), p4)
    draw(PDF(savedir*"p5.pdf"), p5)
    draw(PDF(savedir*"p6.pdf"), p6)
end



# Gadfly.pop_theme()
println("data saved")
