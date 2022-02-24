using JSON
using ColorSchemes
using Compose
import Cairo, Fontconfig
include("../library/abf2_analysis_tool.jl")
include("../library/abf2_pnscan.jl")
println("Compiled.")

fullscan_dir = "/Users/pcs/data/ABF/rawdata/nu2-pn/d1-abf-sf/merged/"
save_dir = "/Users/pcs/data/ABF/analysis/nu2-pn-phase-dis-eff-new/single-summary/"
config_fn = "/Users/pcs/data/ABF/rawdata/nu2-pn/d1-abf-sf/config"

config = JSON.parsefile(config_fn)
p = readconfig(config)
# Read all JLDs in fullscan_dir
arr_full_fn = readdir_jld(fullscan_dir, join = true)
arr_fn = readdir_jld(fullscan_dir, join = false)


E_edges = range(-0.5, 0.5, length = 102)
for i in 1:length(arr_fn)
    p = summary_single_visual(arr_full_fn[i], E_edges,theme)
    println("$(i)/$(length(arr_fn))")
    draw(PNG(save_dir*"$(arr_fn[i]).png", dpi = 150), p)
end
