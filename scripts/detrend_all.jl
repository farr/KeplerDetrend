using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using ArgParse
using DataFrames
using FITSIO
using HDF5
using KeplerDetrend
using Printf
using ProgressLogging
using Statistics

s = ArgParseSettings()
@add_arg_table s begin
    "--cbv_threshold"
        arg_type = Float64
        default = 5.0
        help = "include CBVs until the median basis coefficient significance is below this threshold (in sigma)"
    "--n_cbvs"
        arg_type = Int
        default = 0
        help = "number of CBVs to use for detrending (<= 0 means use var threshold)"
end

parsed_args = parse_args(s)

cbv_threshold = parsed_args["cbv_threshold"] 
n_cbvs = parsed_args["n_cbvs"]

qtrs = 1:17

@progress name="Quarters" for q in qtrs
    basis, svs = h5open(joinpath(@__DIR__, "..", "lc-files", "HLSP", "CBV", @sprintf("Q%02d.h5", q)), "r") do file
        read(file, "basis"), read(file, "singular_values")
    end


    if n_cbvs <= 0 
        fitsdir = joinpath(@__DIR__, "..", "lc-files", "HLSP", @sprintf("Q%02d", q))
        dfs = []
        @progress name="Load Files" for f in readdir(fitsdir)
            if endswith(f, ".fits")
                FITS(joinpath(fitsdir, f), "r") do file
                    df = DataFrame(file[2])
                    df[!, :FLUX_ERR_EST] .= difference_white_noise_estimate(df[!, :FLUX])
    
                    push!(dfs, df)
                end
            end
        end
        nc = num_cbvs_threshold(basis, dfs, threshold=cbv_threshold)
    else
        nc = n_cbvs
    end
    
    M = full_detrend_basis_to_detrend_design_matrix(basis, nc)
    @info "Using $nc CBVs for detrending Q$(q)"

    fitsdir = joinpath(@__DIR__, "..", "lc-files", "HLSP", @sprintf("Q%02d", q))

    @progress name="Files" for f in readdir(fitsdir)
        if endswith(f, ".fits")
            FITS(joinpath(fitsdir, f), "r") do file
                lc = DataFrame(file[2])
                lc[!, :FLUX_ERR_EST] .= difference_white_noise_estimate(lc[!, :FLUX])

                detrend!(lc, M)

                fbase, fext = splitext(f)

                outfile = fbase * "_detrended" * fext

                FITS(joinpath(fitsdir, outfile), "w") do outf
                    header1 = read_header(file[1])
                    header1["NCBV"] = nc
                    set_comment!(header1, "NCBV", "Number of CBVs used to detrend")

                    write(outf, zeros(1,1), header=header1)

                    header2 = read_header(file[2])
                    header2["NAXIS2"] = 8*5
                    header2["TFIELDS"] = 5

                    header2["TTYPE4"] = "FLUX_DETREND"
                    header2["TFORM4"] = "D"
                    header2["TUNIT4"] = "e-/s"
                    set_comment!(header2, "TTYPE4", "Detrended flux")
                    set_comment!(header2, "TUNIT4", "column units: electrons per second")

                    header2["TTYPE5"] = "FLUX_DETREND_ERR_EST"
                    header2["TFORM5"] = "D"
                    header2["TUNIT5"] = "e-/s"
                    set_comment!(header2, "TTYPE5", "Detrended flux error estimate")
                    set_comment!(header2, "TUNIT5", "column units: electrons per second")

                    write(outf, 
                    ["CADENCE", "TIME", "FLUX", "FLUX_DETREND", "FLUX_DETREND_ERR_EST"],
                    [lc[!, :CADENCE], lc[!, :TIME], lc[!, :FLUX], lc[!, :FLUX_DETREND], lc[!, :FLUX_ERR_EST]],
                    header=header2)
                end
            end
        end
    end
end