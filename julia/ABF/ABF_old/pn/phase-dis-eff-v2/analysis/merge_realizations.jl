## Merge dataframes of same file names inside different folders.

include("../library/abf2_analysis_tool.jl")

## set open directories
# Set folders where full scan result for each realizations are stored
root_dir = "/Users/pcs/data/ABF/rawdata/nu2-pn/d1-nu2-sf-size-effect/rawdata/"
load_dir = root_dir.*["R$(i)/" for i in 1:7]

## set save directory
save_dir = root_dir*"merged/"

fns = readdir_jld.(load_dir, join = false)

for i in 1:length(fns[1])
    save_fn = save_dir*fns[1][i]
    dict = JLD.load(load_dir[1]*fns[1][i])
    count = 0
    for j in 2:length(fns)
        idx = findfirst(s-> s == fns[1][i], fns[j])
        if !(isnothing(idx))
            dict_temp = JLD.load(load_dir[j]*fns[j][idx])
            append!(dict["data"], dict_temp["data"])
            count += 1
        else
            println("Warning: a file is missing")
        end
    end
    dict = JLD.save(save_fn, dict)
    println("$(fns[1][i]):  $(count) files merged.")
end
