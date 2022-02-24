#
#  scan.jl
#
#  Copyright 2020 김영준 <yeongjun@ust.ac.kr>
#                 Alexei Andreanov <aalexei@ibs.re.kr>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
import JSON
using MAT
using ArgParse
using Distributed
include("xi_abf2_dist.jl")

function main(args)
#   Default values
#   Command line arguments
    opts = ArgParseSettings(description="Scan and compute localisation lengths for all parameters of nu=2 ABF")
    @add_arg_table! opts begin
    "c"
        help = "configuration"
        #nargs = 1
        arg_type = AbstractString
     "--nprocs"
         help = "Number of cores"
         arg_type = Int
         default = 0
    end

#   Parse the arguments
    popts   = parse_args(opts)
    t = addprocs(popts["nprocs"])
    @everywhere include("xi_abf2_dist.jl")

    config  = JSON.parsefile(popts["c"])
    params = config_read(config)   # create the Parameters struct from Dict
    @time ξ_all = ξ_scan_pmap(params)                 # compute all the localisation lengths

    config["xi_all"]  = ξ_all                                           #
    fn_all = splitext(basename(popts["c"]))                             # output filename
    fn_out = fn_all[1] * "-xi_all.mat"
    matwrite(fn_out, config, compress=true)
    rmprocs(t)
                          # write compressed MATLAB file
end
@time main(ARGS)
