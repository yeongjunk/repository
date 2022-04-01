using ArgParse, JSON
using LinearAlgebra
using MAT
using Plots
inspectdr()
ENV["JULIA_COPY_STACKS"] = 1
include("./diag.jl") # read parameters from configuration file
function main(ARGS)
    opts = ArgParseSettings(description="Scan and compute pn for all parameters of nu=2 ABF")
    @add_arg_table! opts begin
    "c"
        help = "configuration"
        arg_type = AbstractString
    end
    # Parse the arguments
    popts   = parse_args(opts)
    config  = JSON.parsefile(popts["c"])
    p = readconfig(config)
    if p.num_blas != 0
        LinearAlgebra.BLAS.set_num_threads(p.num_blas)
    end
    @time dct = abf_scan(p)
    matwrite("output.mat", dct)
    logpsi = dct["logpsi_mean"]
    logpsi_ste = dct["logpsi_ste"]
    log10_logpsi_ste = dct["log10_logpsi_ste"]
    p1 = plot(logpsi, yerror = logpsi_ste);
    xlabel!("sites")
    ylabel!("g")
    p2 = plot(logpsi.^2, yerror = 2logpsi_ste);
    xlabel!("sites")
    ylabel!("g.^2	")
    p3 = plot(log10.(1:p.L÷4), log10.(logpsi), yerror =log10_logpsi_ste);
    xlabel!("log10(sites)")
    ylabel!("log10(g)")
    savefig(p1, "p1.pdf")  
    savefig(p2, "p2.pdf")
    savefig(p3, "p3.pdf")
end

main(ARGS)
