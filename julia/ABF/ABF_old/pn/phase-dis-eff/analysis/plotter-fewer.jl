using Gadfly
import Cairo, Fontconfig
using CSV
using DataFrames
using ColorSchemes
using LinearAlgebra
using Compose
topdir = "/Users/pcs/data/ABF/rawdata/nu2-pn/phase-dis-eff2/"
topdir2 = "/Users/pcs/data/ABF/rawdata/nu2-pn/phase-dis-eff2/del0/"

savedir = "/Users/pcs/data/ABF/analysis/nu2-pn-phase-dis-eff2/"
fn = "pn_dataframe.csv"
#------- Load data and merge width δ=0 -------#
df = DataFrame(CSV.File(topdir*fn))
df2 = DataFrame(CSV.File(topdir2*fn))

df.E = round.(df.E, digits = 12)
df2.E = round.(df2.E, digits = 12)
rename!(df2, Dict(:Delta => "L10_Delta"))

df.L10_Delta = string.(df.L10_Delta)
df2.L10_Delta = string.(df2.L10_Delta)
df2.L10_Delta = replace.(df2.L10_Delta, "0.0" => "δ = 0")
df = vcat(df2,df)

#------- Load data and merge width δ=0 -------#
dfe0 = df[df.E .==0.,:]
dfe0_2 = dfe0[(dfe0.L10_Delta .=="-2.0") .| (dfe0.L10_Delta .=="-1.0").| (dfe0.L10_Delta .=="0.0").| (dfe0.L10_Delta .=="1.0") .| (dfe0.L10_Delta .=="1.0") .|  (dfe0.L10_Delta .=="2.0") .| (dfe0.L10_Delta .=="3.0") , :]

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
    panel_line_width = 2pt ,
    plot_padding = [1mm, 5mm, 1mm, 1mm],
    grid_line_width = 1pt,
    grid_line_style = :solid,
)

# Gadfly.with_theme(theme) do
# #------------------------------ PN vs THETA ------------------------------#
#     set_default_plot_size(11cm, 9cm)
#     xticks = collect(0.05:0.05:0.25)
#     pushfirst!(xticks, 0.01)
#     yticks = 1:2:11
#     l1 = layer(dfe0[dfe0.L10_Delta .== "δ = 0",:], x = :Theta, y = :PN_mean, Geom.line, color =[colorant"black"],
#         style(line_width = 0.5mm, point_size = 2.5pt))
#     l2 = layer(dfe0[dfe0.L10_Delta .== "-1.0",:], x = :Theta, y = :PN_mean, Geom.line, color =[colorant"purple"],
#     style(line_width = 0.5mm, point_size = 2.5pt))
#
#     l3 = layer(dfe0[dfe0.L10_Delta .== "0.0",:], x = :Theta, y = :PN_mean, Geom.line, color =[colorant"blue"],
#     style(line_width = 0.6mm, point_size = 3pt))
#     l4 = layer(dfe0[dfe0.L10_Delta .== "4.0",:], x = :Theta, y = :PN_mean, Geom.line, color =[colorant"red"],
#     style(line_width = 0.5mm, point_size = 2.5pt))
#     l5 = layer(dfe0[dfe0.L10_Delta .== "5.0",:], x = :Theta, y = :PN_mean, Geom.line, color =[colorant"orange"],
#     style(line_width = 0.5mm, point_size = 2.5pt))
#     p1 = plot(l1,l2,l3,l4,l5,  style(highlight_width=0mm),
#     Guide.xticks(ticks=xticks),
#     Guide.yticks(ticks=yticks),
#     Guide.xlabel("Θ/π"),
#     Guide.ylabel("PN")
#     )
#     display(p1)
# end


df = DataFrame(CSV.File(topdir*fn))
df.E = round.(df.E, digits = 12)
dfe0 = df[df.E .==0.,:]
Gadfly.with_theme(theme) do
    xticks = collect(0.00:0.05:0.25)
    yticks = 1:2:11
    p1 = plot(dfe0[dfe0.L10_Delta .== 0.0,:], x =:Theta, y =:PN_mean, color=[colorant"red"], Geom.line,
        Guide.yticks(ticks=yticks),
        Guide.xticks(ticks=xticks),
        Guide.xlabel("Θ/π"),
        Guide.ylabel("PN"),
        style(line_width=0.6mm),
        Guide.annotation(compose(context(), text(0.3w, 0.7h, "δ = 1.0"), fontsize(18pt)))
        )
        display(p1)

#------------------------------ PN vs delta at Θ=0.25π ------------------------------#
    xticks = -3:1:5
    yticks = 1:2:11
    p2 = plot(dfe0[dfe0.Theta .== 0.25 ,:], x =:L10_Delta, y =:PN_mean, color=:Theta, Geom.line,
        Guide.yticks(ticks=yticks),
        Guide.xticks(ticks=xticks),
        Guide.xlabel("log<sub>10</sub>[δ]"),
        Guide.ylabel("PN"),
        Guide.colorkey(title = "Θ/π",pos=[8cm,-9cm]),
        style(line_width=0.6mm),
        Scale.color_discrete_manual("red","purple","green"),
        )
    display(p2)
#------------------------------ SAVE------------------------------#
    draw(PDF(savedir*"p1.pdf"), p1)
    draw(PDF(savedir*"p2.pdf"), p2)
end
