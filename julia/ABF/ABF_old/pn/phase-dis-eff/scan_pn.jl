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
    end

    # Parse the arguments
    popts   = parse_args(opts)
    config  = JSON.parsefile(popts["c"])

    params = readconfig(config)   # create the Parameters struct from Dict config
    if isclean(params) && p.E_abs_bw == false
        error("Clean ABF is included which has zero bandwidth. You cannot compute relative bandwidth.")
    end

    if params.W == 0 && (params.core == 1 || params.core == 2)
        error("For the scale free model, W cannot be zero.")
    end
    
    BLAS.set_num_threads(1) #Turn off BLAS multi-threading
    time_begin = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    # @time N, PN_mean, PN_var, PN_hist, E_bw = scan_pn(params) # Two lines for test
    if params.core == 1
        f = pn_sf_ph_on
    elseif params.core == 2
        f = pn_sf_ph
    elseif params.core == 3
        f = pn_fd_ph_on
    elseif params.core == 4
        f = pn_fd_ph
    end

    @time N, PN_mean, PN_var, PN_hist, E_bw = scan_pn(f,params)

    time_end = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    info = "Beign:"*time_begin*"\n"*"End: "*time_end*"\n"

    config["num_samp"] = N
    config["PN_mean"]  = PN_mean
    config["PN_var"]  = PN_var
    config["PN_hist"]  = PN_hist
    config["E_bw"] = E_bw
    config["info"] = info

    fn_out = "pn.jld"
    save(fn_out, config)
    println("DATA SAVED")
    println("max(PN_mean) = $(maximum(PN_mean))")

end

main(ARGS)
