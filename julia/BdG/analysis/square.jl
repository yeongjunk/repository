using CSV
using Plots
using DataFrames
using LaTeXStrings

function relative(f, r; sp)
    p = plot!()
    lims = f(p[sp])
    return lims[1] + r * (lims[2]-lims[1])
end
relativex(r; sp::Int=1) = relative(Plots.xlims, r; sp=sp)
relativey(r; sp::Int=1) = relative(Plots.ylims, r; sp=sp)


dir = "/Users/kimyeongjun/data/BdG/square/"

fn_2 = ["bdgsquare-W4.0-ga2.01-nev1-R50-L40.csv"
    "bdgsquare-W4.0-ga2.01-nev1-R50-L80.csv"
    "bdgsquare-W4.0-ga2.01-nev1-R50-L120.csv"
    "bdgsquare-W4.0-ga2.01-nev1-R50-L160.csv"
    "bdgsquare-W4.0-ga2.01-nev1-R50-L200.csv"]

fn_5=["bdgsquare-W4.0-ga5.0-nev1-R50-L40.csv"
    "bdgsquare-W4.0-ga5.0-nev1-R50-L80.csv"
    "bdgsquare-W4.0-ga5.0-nev1-R50-L120.csv"
    "bdgsquare-W4.0-ga5.0-nev1-R50-L160.csv"
    "bdgsquare-W4.0-ga5.0-nev1-R50-L200.csv"
    "bdgsquare-W4.0-ga5.0-nev1-R50-L500.csv"]


fn_20 = ["bdgsquare-W4.0-ga20.0-nev1-R50-L40.csv"
    "bdgsquare-W4.0-ga20.0-nev1-R50-L80.csv"
    "bdgsquare-W4.0-ga20.0-nev1-R50-L120.csv"
    "bdgsquare-W4.0-ga20.0-nev1-R50-L160.csv"
    "bdgsquare-W4.0-ga20.0-nev1-R50-L200.csv"]

fn_40 = [
    "bdgsquare-W4.0-ga40.0-nev1-R50-L40.csv"
    "bdgsquare-W4.0-ga40.0-nev1-R50-L80.csv"
    "bdgsquare-W4.0-ga40.0-nev1-R50-L120.csv"
    "bdgsquare-W4.0-ga40.0-nev1-R50-L160.csv"
    "bdgsquare-W4.0-ga40.0-nev1-R50-L200.csv"]

fn_2 = dir.*fn_2
fn_5 = dir.*fn_5
fn_20 = dir.*fn_20
fn_40 = dir.*fn_40

df_2 = [CSV.read(fn_2[i], DataFrame) for i in 1:length(fn_2)]
df_5 = [CSV.read(fn_5[i], DataFrame) for i in 1:length(fn_5)]
df_20 = [CSV.read(fn_20[i], DataFrame) for i in 1:length(fn_20)]
df_40 = [CSV.read(fn_40[i], DataFrame) for i in 1:length(fn_40)]

default(
    framestyle = :box,
    size = (600,400),
    # right_margin = [3mm 0mm],
    grid = false,
    minorticks = true,
    legend = :bottomleft,
    fontfamily = "computer modern",
    tickfontsize = 12,
    guidefontsize = 12,
    legendfontsize = 10,
    annotationfontsize = 12, palette = :default)
marker = (:circle, 3., 1., stroke(-0.5, 1., :black))
line = (:line, :solid, 2.5)

L = ["40" "80" "120" "160" "200" "500"]
L_int = [40, 80, 120, 160, 200, 500]
function plot_pn(df, n)
    p1 = plot(df[1].λ, df[1].pn, yaxis = :log10, m = marker, line = line, label = L"L = %$(L[1])")
    for i in 2:length(df)
        plot!(p1, df[i].λ, df[i].pn, yaxis = :log10, m = marker, line = line, label = L"L = %$(L[i])")
    end
    vline!(p1, [sqrt(8*n)], label = L"\sqrt{8ga}", c = :black)
    xlabel!(L"\lambda")
    ylabel!(L"P")
    annotate!(p1, sp=1,[(relativex(0.5; sp=1), relativey(0.6; sp=1), L"ga = %$n", 25)])
    return p1
end



p = plot_pn(df_2, 2)
p1 = plot_pn(df_5, 5)
p2 = plot_pn(df_20, 20)
p3 = plot_pn(df_40, 40)
p_full = plot(p,p1,p2,p3, size = (1200, 800))
savefig(p_full, dir*"figures/ga_full.pdf")
log_L = log.(L_int)
df_full = [df_2, df_5, df_20, df_40]
n = [2, 5, 20, 40]
for j in 1:length(df_full) 
df_n = df_full[j]
tau = Array{Float64}[]
for i in 1:length(df_n)
    push!(tau, log.(df_n[i].pn)./log_L[i])
end    

p_tau = plot()
for i in 1:length(tau)
    plot!(p_tau, df_n[i].λ, tau[i], line = line, marker = marker, label = "L = $(L_int[i])")
end
ylabel!(L"\tau")
annotate!(p_tau, sp=1, [(relativex(0.5;sp=1), relativey(0.6;sp=1), L"ga = %$(n[j])", 25)])
savefig(p_tau, dir*"tau$(j).pdf")
end

