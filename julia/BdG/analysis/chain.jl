using CSV
using Plots
using DataFrames
using LaTeXStrings
dir = "/Users/pcs/data/BdG/chain/"
fn_2 = dir*"bdgchain-W4.0-ga2.01-nev1-R50.csv"
fn_3 = dir*"bdgchain-W4.0-ga3.0-nev1-R50.csv"
fn_5 = dir*"bdgchain-W4.0-ga5.0-nev1-R50.csv"



df_2 = CSV.read(fn_2, DataFrame)
df_3 = CSV.read(fn_3, DataFrame)
df_5 = CSV.read(fn_5, DataFrame)

default(
    framestyle = :box,
    size = (600,400),
    # right_margin = [3mm 0mm],
    grid = false,
    minorticks = true,
    legend = :topright,
    fontfamily = "computer modern",
    tickfontsize = 12,
    guidefontsize = 12,
    legendfontsize = 12,
    annotationfontsize = 12, palette = :default)
marker = (:circle, 4., 1., stroke(-0.5, 1., :black))
line = (:line, :solid, 2.5)

p = plot(df_2.λ, df_2.pn, yaxis = :log10, m = marker, line = line, label = "ga = 2")
plot!(p, df_3.λ, df_3.pn, yaxis = :log10, m = marker, line = line, label = "ga = 3")
plot!(p, df_5.λ, df_5.pn, yaxis = :log10, m = marker, line = line, label = "ga = 5")
xlabel!(L"\lambda")
ylabel!(L"P")

savefig(p, dir*"figure.pdf")
