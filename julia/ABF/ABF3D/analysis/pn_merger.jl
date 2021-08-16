    using CSV, DataFrames
θ = range(0.0001, 0.25, length = 26)

##
L = repeat([30], 10)
R = repeat([5], 10)
dir = ["/Users/pcs/data/ABF3D/full-ipr/L30/$i/" for i in 1:10]



savedir = "/Users/pcs/data/ABF3D/full-ipr/L30/merged/"
##
function load_file(dir, i, j, L, R, θ)
    fn = "L$(L[i])_Th$(j)_R$(R[i]).csv"
    return CSV.read(dir*fn, DataFrame)
end

i = 1
for j in 1:26
    df = [load_file(dir[k], i, j, L[k], R[k], θ) for k in 1:10]

    for k in 2:10
        df[k].r .+= df[k-1][end,2]
    end
    df_new = df[1]
    for k in 2:10
        append!(df_new, df[k])
    end
    if issorted(df_new.r)
        CSV.write(savedir*"L$(L[i])_Th$(j)_R$(df_new[end,2]).csv", df_new)
    end
end
