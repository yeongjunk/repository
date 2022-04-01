using JSON, ArgParse
using CSV
include("./xi.jl")
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "config"
            help = "configuration file"
            required = true
        "--out"
            help = "Output file name"
            arg_type = String
            default = "xi.csv"
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()
    display(args)
    config = JSON.parsefile(args["config"])
    p = read_config(config)
    xi = scan_xi(p)
    CSV.write(args["out"], DataFrame(E = p.E, xi = vec(xi)))
end

main()
