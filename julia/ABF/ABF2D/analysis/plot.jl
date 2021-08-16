using Plots, MAT
using Plots.PlotMeasures
#-------------- ROAG --------------#
L = [20 40 60 80 100]
R = [10000 2500 1000 1000 500]
θ = range(0.0001, 0.25, length = 26)
E_edges = collect(0.65:0.1:1.35)
pushfirst!(E_edges,0.60)
push!(E_edges,  1.4)

r_data = matread(savedir*"roag.mat")
r_m_m = r_data["r_m_m"]
marker = (:circle, 4, 1., stroke(-0.5, 1., :black))
line = (:line, :solid, 2)
palette_roag = :Dark2_5
default(
    framestyle = :box,
    size = (600,400),
    # right_margin = [3mm 0mm],
    grid = false,
    minorticks = true,
    legend = (0.1, 0.75),
    fontfamily = "computer modern",
    tickfontsize = 13,
    guidefontsize = 13,
    legendfontsize = 13, palette = :default)

p = plot(
    xlabel = L"$\Theta / \pi$",
    ylabel = L"$\langle r \rangle $")
for i in 1:5
    plot!(θ, r_m_m[5, :, i], line = line, marker = marker, label = L"$L = %$(L[i])$")
end
annotate!(0.03, 0.52, L"E = 1")
p
hline!([0.3863], c = :red, ls = :dash, label = "Poisson")
hline!([0.5307], c = :blue, ls = :dash, label = "GOE")


p1 = plot(
    xlabel = L"$1/L$",
    ylabel = L"$\langle r \rangle $")
plot!(p1, vec(L).^-1,  vec(r_m_m[5,end,:]), legend = false, line = line, marker = marker, frame = :box)
annotate!(0.015, 0.50, L"\Theta /\pi = 0.25")
## PN
L = [20 40 60]
R = [10000 2500 1000]
θ = range(0.0001, 0.25, length = 26)
E_edges = collect(0.65:0.1:1.35)
pushfirst!(E_edges,0.60)
push!(E_edges,  1.4)

pn_data = matread(savedir*"pn.mat")
pn_m_m = pn_data["pn_m_m"]
pn_exp_diff = pn_data["pn_exp_diff"]
palette_roag = :Dark2_5
default(
    framestyle = :box,
    size = (600,400),
    # right_margin = [3mm 0mm],
    grid = false,
    minorticks = true,
    legend = (0.1, 0.75),
    fontfamily = "computer modern",
    tickfontsize = 13,
    guidefontsize = 13,
    legendfontsize = 13, palette = :tab10)

p2 = plot(
    xlabel = L"$\Theta / \pi$",
    ylabel = L"$\langle PN \rangle $")
for i in 1:3
    plot!(θ, pn_m_m[5, :, i], line = line, marker = marker, label = L"$L = %$(L[i])$")
end
annotate!(0.02, 90, L"E = 1")

p3 = plot(
    xlabel = L"$\Theta / \pi$",
    ylabel = L"$\alpha$")
for i in 1:2
    plot!(p3, θ,  pn_exp_diff[5,:,i], label = L"$L = %$(L[i])$" , line = line, marker = marker, frame = :box)
end

p_full = plot(p, p1, p2, p3, right_margin = 3mm, left_margin = 3mm, layout = 4, size = (1200, 800))

savefig(p_full, savedir*"2D-results.pdf")
