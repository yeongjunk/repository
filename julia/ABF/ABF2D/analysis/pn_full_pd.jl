using DataFrames, CSV, MAT
using ROAG, Binning
using Statistics
using StatsBase
using Interpolations, LaTeXStrings
using Plots
using Glob

##
L = collect(20:20:140)
θ = range(0.0001, 0.25, length = 26)

##
dir = "/Users/pcs/data/ABF2D/pn-pd-sf/"

dirs = [dir*"L20"; dir*"L40"; dir*"L60"; dir*"L80"; dir*"L100"; dir*"L120"; dir*"L140"]

#
##
df = [DataFrame(E = Float64[], PN = Float64[]) for i in 1:length(dirs)]
for j in 1:length(dirs)
    files = glob("*.csv", dirs[j])
    for i in 1:length(files)
        df_temp = CSV.read(files[i], DataFrame)
        append!(df[j], df_temp)
    end
end

del_E = 0.02
e_c = 1.
e_min = e_c - del_E
e_max = e_c + del_E
function filt_E(x, e_min, e_max)
     return (x < e_max) && (x > e_min)
 end


pn_mean = Float64[];
for i in 1:length(df)
    push!(pn_mean, mean(df[i][filt_E.(df[i].E, e_min, e_max), :PN]))
end
log_pn = log.(pn_mean)
dlog_dL = diff(log_pn)./diff(log.(L))
scatter(log.(L), dlog_dL)
