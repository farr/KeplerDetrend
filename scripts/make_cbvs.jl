using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using DataFrames
using HDF5
using FITSIO
using KeplerDetrend
using Printf
using ProgressMeter

qtrs = 1:17

@showprogress for q in qtrs
    fitsdir = joinpath(@__DIR__, "..", "lc-files", "HLSP", @sprintf("Q%02d", q))
    dfs = []
    for f in readdir(fitsdir)
        if endswith(f, ".fits")
            FITS(joinpath(fitsdir, f), "r") do file
                lc = DataFrame(file[2])
                lc[!, :FLUX_ERR_EST] .= difference_white_noise_estimate(lc[!, :FLUX])
                push!(dfs, lc)
            end
        end
    end

    U, S = full_detrend_basis_and_singular_values(dfs)
    
    h5open(joinpath(@__DIR__, "..", "lc-files", "HLSP", "CBV", @sprintf("Q%02d.h5", q)), "w") do file
        file["basis", compress=3, shuffle=()] = U
        file["singular_values", compress=3, shuffle=()] = S
    end
end