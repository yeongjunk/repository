using Plots, MAT
using Plots.PlotMeasures
using LaTeXStrings
using StatsBase

savedir = "/Users/pcs/data/ABF2D/onsite/analyzed-data/"
savedir1d = "/Users/pcs/data/ABF1D/onsite/analyzed-data/"

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

## 1D
L = [100 1000]
R = [10000 1000]
θ = range(0.0001, 0.25, length = 26)
E_edges = collect(0.65:0.1:1.35)
pushfirst!(E_edges,0.60)
push!(E_edges,  1.4)

pn_data = matread(savedir1d*"pn.mat")
pn_m_m = pn_data["pn_m_m"]

p1 = plot(
    xlabel = L"$\theta / \pi$",
    ylabel = L"$\langle PN \rangle $",
    xlims = (0, 0.25))
for i in 1:2
    plot!(p1, θ, pn_m_m[5, :, i], line = line, marker = marker, label = L"$L = %$(L[i])$")
end
p1

## p1_inset
plot!(p1, sp = 2, palette = :default, legend = false,
    inset = (1, bbox(0.04, 0.25, 0.35, 0.3, :bottom, :right)),
    xlabel =  L"$E$",
    ylabel = L"$\langle PN \rangle$",
    xlims = (0.6, 1.4),
    ylims = (2.0, 3.0),
    xticks = [0.6, 1.0 ,1.4],
    yticks = [2.0, 2.5, 3.0])

for i in 1:2
    plot!(p1, sp = 2, midpoints(E_edges), pn_m_m[:, end, i], line = line, marker = marker, label = L"$L = %$(L[i])$")
end
display(p1)
## 2D PN
## PN
L = [20 40 60 80 100]
R = [10000 2500 1000 100 40]
θ = range(0.0001, 0.25, length = 26)
E_edges = collect(0.65:0.1:1.35)
pushfirst!(E_edges,0.60)
push!(E_edges,  1.4)

pn_data = matread(savedir*"pn.mat")
pn_m_m = pn_data["pn_m_m"]
pn_exp_diff = pn_data["pn_exp_diff"]
palette_pn = :Greens_5

p2 = plot(
    palette = palette_pn,
    xlabel = L"$\theta / \pi$",
    ylabel = L"$\langle PN \rangle $",
    xlims = (0, 0.25))
for i in 1:length(L)
    plot!(θ, pn_m_m[5, :, i], line = line, marker = marker, label = L"$L = %$(L[i])$")
end

p_diff = plot(
    xlabel = L"$\theta / \pi$",
    ylabel = L"$\alpha$")
for i in 1:2
    plot!(p_diff, θ,  pn_exp_diff[5,:,i], label = L"$L = %$(L[i])$" , line = line, marker = marker, frame = :box)
end

##
    pn = pn_m_m[5,:,:]
    p3 = plot(
        xlabel = L"$\ln{L}$",
        ylabel = L"$\ln{\langle PN \rangle} $",
        xlims = (log(L[1]), log(L[end])+0.0001),
        xticks = (log.(vec(L))[1:end-1], "ln".*string.(vec(L))[1:end-1]))
        plot!(log.(vec(L)), log.(pn'), line = line,marker = marker, label = :none, line_z = θ',marker_z = θ', c = :viridis, colorbar = false)
    # annotate!(3.27, 4.4, L"$\theta/\pi = 0.25$")
    # annotate!((log(L[1]) + log(L[2]))/2, (log(pn[1]) + log(pn[2]))/2, L"\alpha_{20 \rightarrow 40} ")
    # annotate!((log(L[2]) + log(L[3]))/2, (log(pn[2]) + log(pn[3]))/2, L"\alpha_{40 \rightarrow 60} ")
    # annotate!((log(L[3]) + log(L[4]))/2, (log(pn[3]) + log(pn[4]))/2, L"\alpha_{60 \rightarrow 80} ")
    # annotate!((log(L[4]) + log(L[5]))/2, (log(pn[4]) + log(pn[5]))/2, L"\alpha_{80 \rightarrow 100} ")


    p3

#-------------- ROAG --------------#

L = [20 40 60 80 100]
R = [10000 2500 1000 1000 500]
θ = range(0.0001, 0.25, length = 26)
E_edges = collect(0.65:0.1:1.35)
pushfirst!(E_edges,0.60)
push!(E_edges,  1.4)

r_data = matread(savedir*"roag.mat")
r_m_m = r_data["r_m_m"]
palette_roag = :PuBu_5

p4 = plot(
    palette = palette_roag,
    xlabel = L"$\theta / \pi$",
    ylabel = L"$\langle r \rangle $",
    yticks = 0.40:0.04:0.52,
    xlims = (0, 0.25))
for i in 1:5
    plot!(θ, r_m_m[5, :, i], line = line, marker = marker, label = L"$L = %$(L[i])$")
end
# annotate!(0.12, 0.46, L"\theta /\pi = 0.25")

p4
hline!([0.3863], c = :blue, ls = :dash, label = :none, lw = 2)

p4_1 = plot(
    xlabel = L"$1/L$",
    ylabel = L"$\langle r \rangle $")
plot!(p4_1, vec(L).^-1,  vec(r_m_m[5,end,:]), legend = false, line = line, marker = marker, frame = :box)
annotate!(0.015, 0.52, L"\theta /\pi = 0.25", 25)

##
p_full = plot(p1, p2, p4, p3, layout = 4, size = (1200, 800))
plot!(p3, left_margin = 5mm, right_margin = 7mm)
plot!(p4, left_margin = 5mm, bottom_margin = 3mm)
plot!(p2, left_margin = 2mm)
p_full = plot(p1, p2, p4, p3, layout = 4 ,size = (1200, 800))

annotate!(p_full, sp=1,[(relativex(0.03; sp=1), relativey(0.93; sp=1), text("(a)",:black, :left, 22, "computer modern"))])
annotate!(p_full, sp=3,[(relativex(0.03; sp=3), relativey(0.93; sp=3), text("(b)",:black ,:left, 22, "computer modern"))])
annotate!(p_full, sp=4,[(relativex(0.03; sp=4), relativey(0.93; sp=4), text("(c)",:black ,:left, 22, "computer modern"))])
annotate!(p_full, sp=5,[(relativex(0.03; sp=5), relativey(0.93; sp=5), text("(d)",:black ,:left, 22, "computer modern"))])

annotate!(p_full, sp=1,[(relativex(0.5; sp=1), relativey(0.93; sp=1), L"\mathbf{\textit{d=1}}", 22)])
annotate!(p_full, sp=3,[(relativex(0.5; sp=3), relativey(0.93; sp=3), L"\mathbf{\textit{d=2}}", 22)])
annotate!(p_full, sp=4,[(relativex(0.5; sp=4), relativey(0.93; sp=4), L"\mathbf{\textit{d=2}}", 22)])
annotate!(p_full, sp=5,[(relativex(0.5; sp=5), relativey(0.93; sp=5), L"\mathbf{\textit{d=2}}", 22)])


savefig(p_full,savedir1d*"full_1d2d_2.pdf")

p_full
