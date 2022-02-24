using JLD
using JSON
using ArgParse
using Dates
include("library/abf2_pnscan.jl")

function main(args)
    opts = ArgParseSettings(description="Scan and compute pn for all parameters of nu=2 ABF")
    @add_arg_table! opts begin
    "c"
        help = "configuration"
        arg_type = AbstractString

    "--log"
        help = "make scale W log i.e. 10^W"
        action = :store_true
    end

    # Parse the arguments
    popts   = parse_args(opts)
    config  = JSON.parsefile(popts["c"])
    set_logscale  = popts["log"]

    params = readconfig(config)   # create the Parameters struct from Dict config
    BLAS.set_num_threads(1) #Turn off BLAS multi-threading
    time_begin = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    @time smp_num, pn_mean, pn_var, pn_hist = scan_pn(params; logscale = set_logscale)
    time_end = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    info = "Beign:"*time_begin*"\n"*"End: "*time_end*"\n"*"logscale: $(set_logscale)"

    config["smp_num"] = smp_num
    config["pn_mean"]  = pn_mean
    config["pn_var"]  = pn_var
    config["pn_hist"]  = pn_hist
    config["info"] = info

    fn_out = "pn.jld"
    save(fn_out, config) # write compressed MATLAB file
end

main(ARGS)
