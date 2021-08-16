using Gadfly; import Cairo, Fontconfig
using CSV
using DataFrames
using ColorSchemes
using LinearAlgebra
using Compose
topdir = "/Users/pcs/data/ABF/rawdata/nu2-pn/phase-dis-eff2/hist-csv/"
savedir = "/Users/pcs/data/ABF/analysis/nu2-pn-phase-dis-eff2/"

#-----------------------------PLOT HISTOGRAMS-----------------------------#
for j in 1:4:25
    fn = "th$(j)-del1.csv"
    df = DataFrame(CSV.File(topdir*fn))
    df = df[df.E .== 0.,:]
    df.Num = normalize(df.Num)

    for i in 11:10:70
        fn = "th$(j)-del$(i).csv"
        df_temp = DataFrame(CSV.File(topdir*fn))
        df_temp = df_temp[df_temp.E .== 0.,:]
        df_temp.Num = normalize(df_temp.Num)
        df = vcat(df, df_temp)
    end
    df.L10_Delta = string.(log10.(df.L10_Delta))
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
        panel_stroke    =   "black",
        panel_line_width = 2pt ,
        plot_padding = [1mm, 5mm, 1mm, 1mm],
        line_width = 1pt,
        grid_line_width = 1pt,
        grid_line_style = :solid,
    )
    theme_big = Theme(
        panel_fill= "white",
        #-----FONT TYPE-----#
        major_label_font= ft,
        minor_label_font= ft,
        key_title_font  = ft,
        key_label_font  = ft,
        #-----FONT SIZE-----#
        major_label_font_size=10pt,
        minor_label_font_size=10pt,
        key_title_font_size=9pt,
        key_label_font_size=9pt,
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
    set_default_plot_size(13cm, 10cm)

    xticks = [1, 2, 3, 5, 10, 15, 20]
    Gadfly.with_theme(theme) do
        p1 = plot(df, x = :PN, y = :Num, color = :L10_Delta, Geom.line,
            Guide.xlabel("PN"),
            Guide.ylabel("Probability Amplitude", orientation = :vertical),
            Guide.colorkey(title = "log<sub>10</sub>[δ]",pos=[8cm,-1.1cm]),
            Guide.xticks(ticks = xticks),
            style(line_width = 0.6mm),
            Guide.annotation(compose(context(), text(0.3w, 0.2h, "Θ = $(df.Theta[1])π"), fontsize(18pt)))
        );
        display(p1)

        p2 = plot(df, x = :PN, y = :Num, color = :L10_Delta, Geom.line,
            Guide.xlabel("PN"),
            Guide.ylabel("Probability Amplitude", orientation = :vertical),
            Guide.colorkey(title = "log<sub>10</sub>[δ]"),
            Guide.xticks(ticks = xticks),
            style(line_width = 0.6mm),
            Scale.y_log10,
            Guide.annotation(compose(context(), text(0.3w, 0.2h, "Θ = $(df.Theta[1])π"), fontsize(18pt)))
        );

        draw(PDF(savedir*"hist-th$(j).pdf"), p1)
        draw(PDF(savedir*"hist-th$(j)-log.pdf"), p2)
    end
end

#-----------------------------PLOT BANDWIDTH-----------------------------#
topdir = "/Users/pcs/data/ABF/rawdata/nu2-pn/phase-dis-eff2/"
fn2 = "pn_bw.csv"
df = DataFrame(CSV.File(topdir*fn2))
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
    xticks = 0.00:0.05:0.25
    p1 = plot(df, x = :Theta, y = :Δ, color=:L10_Delta, Geom.line,
        Guide.xticks(ticks = xticks),
        Scale.ContinuousColorScale(p -> get(ColorSchemes.jet, p)),
        Guide.xlabel("Θ/π"),
        Guide.ylabel("Δ"),
        Guide.colorkey(title="log<sub>10</sub>[δ]"),
        style(line_width=0.6mm)
    );
    display(p1)

    xticks = 0.01:0.01:0.05
    p2 = plot(df[df.Theta .<= 0.05,:], x = :Theta, y = :Δ, color=:L10_Delta, Geom.line,
        Guide.xticks(ticks = xticks),
        Coord.cartesian(ymin=0, ymax=2),
        Scale.ContinuousColorScale(p -> get(ColorSchemes.jet, p)),
        Guide.xlabel("Θ/π"),
        Guide.ylabel("Δ"),
        Guide.colorkey(title="log<sub>10</sub>[δ]"),
        style(line_width=0.6mm)
    );
    display(p2)

    xticks = -2:1:5
    p3 = plot(df, x = :L10_Delta, y = :Δ, color=:Theta, Geom.line,
        Guide.xticks(ticks = xticks),
        Scale.ContinuousColorScale(p -> get(ColorSchemes.jet, p)),
        Guide.xlabel("log<sub>10</sub>[δ]"),
        Guide.ylabel("Δ"),
        Guide.colorkey(title="Θ/π"),
        style(line_width=0.6mm),
        Scale.y_log10
    );
    display(p3)

    draw(PDF(savedir*"bw1.pdf"), p1)
    draw(PDF(savedir*"bw2.pdf"), p2)
    draw(PDF(savedir*"bw3.pdf"), p3)
end
