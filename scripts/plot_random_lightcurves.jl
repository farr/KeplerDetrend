using Pkg
Pkg.activate("..")

using CairoMakie
using DataFrames
using FITSIO
using LaTeXStrings
using Printf
using Statistics

function plot_random_lightcurve(; qtrs = 1:17, rescale_lc = true)
    q = rand(qtrs)
    fitsdir = joinpath(@__DIR__, "..", "lc-files", "HLSP", @sprintf("Q%02d", q))
    files = readdir(fitsdir)
    detrend_files = filter(f -> endswith(f, "_detrended.fits"), files)
    f = rand(detrend_files)

    lc = DataFrame(FITS(joinpath(fitsdir, f), "r")[2])

    if rescale_lc
        flux_mean = mean(lc[!, :FLUX])
        flux = lc[!, :FLUX_DETREND] ./ flux_mean
        flux_err = lc[!, :FLUX_DETREND_ERR_EST] ./ flux_mean
    else
        flux = lc[!, :FLUX_DETREND]
        flux_err = lc[!, :FLUX_DETREND_ERR_EST]
    end
    fig = Figure()
    a = Axis(fig[1, 1], title=f, xlabel=L"t / \mathrm{d}", ylabel=L"f / \mathrm{e}^- \, \mathrm{s}^{-1}")
    errorbars!(a, lc[!, :TIME], flux, flux_err)
    fig
end