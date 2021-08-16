using JSON
using ArgParse

@doc """
Create configurations file according to specified indices
"""
function config_div(configfile::String, idx::Int64...)
    config  = JSON.parsefile(configfile) #open JSON file(as dictionary)
    mkdir("Config_div") #divided configuration files saved here

    n = length(idx)
    @assert idx[n] == config["W_num"] #conditions that integers should satisfy
    @assert idx[1] == 1

    for i in 1:n-1
        if i > 1
            @assert idx[i] > idx[i-1]
        end
        if i == n-1
            config["W_first"] = idx[i]
            config["W_last"] = idx[i+1]
        else
            config["W_first"] = idx[i]
            config["W_last"] = idx[i+1] - 1
        end
        open("Config_div/"*configfile*"_$(i)", "w") do f # save configuration files
            JSON.print(f, config, 4)
        end
    end
end

function main(args)
#   Default values
#   Command line arguments
    opts = ArgParseSettings(description="Scan and compute localisation lengths for all parameters of nu=2 ABF")
    @add_arg_table! opts begin
    "c"
        help = "configuration"
        #nargs = 1
        arg_type = AbstractString
    end

#   Parse the arguments
    popts   = parse_args(opts)
    config_div(popts["c"],1,11,21,31,41,50)
end

main(ARGS)
