using CSV, DataFrames
θ = range(0.0001, 0.25, length = 26)

##
L_1 = 10
R_1 = 2000
dir_1 = "/Users/pcs/data/ABF3D/full-e/L10/"

L_2 = 10
R_2 = 4000
dir_2 = "/Users/pcs/data/ABF3D/full-e/L10/2/"


savedir = "/Users/pcs/data/ABF3D/full-e/L10/merged/"
##
function load_file(dir, i, j, L, R, θ)
    fn = "L$(L[i])_Th$(j)_R$(R[i]).csv"
    return CSV.read(dir*fn, DataFrame)
end

i = 1
for j in 1:26
    df_1 = load_file(dir_1, i, j, L_1, R_1, θ)
    df_2 = load_file(dir_2, i, j, L_2, R_2, θ)

    println(df_1[end,2])
    df_2.r .+= df_1[end,2]
    println(df_2[end,2])
    append!(df_1, df_2)

    CSV.write(savedir*"L$(L_1[i])_Th$(j)_R$(df_1[end,2]).csv", df_1)
end
