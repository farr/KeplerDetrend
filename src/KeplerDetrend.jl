module KeplerDetrend

using DataFrames
using LinearAlgebra
using Statistics

export difference_white_noise_estimate, full_detrend_basis_to_detrend_design_matrix, full_detrend_basis_and_singular_values, detrend!, num_cbvs_chi2_threshold

@doc raw"""
    difference_white_noise_estimate(flux)

Returns an estimate of the white noise level in the flux time series using a
method that is insensitive to long-term trends.

The estimate is computed as 
```math
\sqrt{\frac{1}{2} \text{var}(\Delta f)}
```
"""
function difference_white_noise_estimate(flux)
    df = diff(flux)
    sqrt(var(df)/2)
end

@doc raw"""
    full_detrend_basis_and_singular_values(data_frames; flux_col_name="FLUX", flux_err_col_name="FLUX_ERR_EST")

Returns a detrending basis (SVD U-matrix) and singular values (S-matrix) for the
set of lightcurves in the given data frames.

Gives a set of basis vectors (i.e. the full U-matrix of the SVD) and singular
values (i.e. the S-matrix of the SVD) for the set of lightcurves in the given
data frames.  Each data frame should have a column for the flux and a colum for
the (estimated) flux uncertainty; by default these columns are named `FLUX` and
`FLUX_ERR_EST`.  

The detrending basis is constructed by SVD of the matrix whose columns are the
zero-mean, whitened lightcurves.  This implies that the constant vector is in
the nullspace of the SVD.

"""
function full_detrend_basis_and_singular_values(dfs; flux_col_name="FLUX", flux_err_col_name="FLUX_ERR_EST")
    M = hcat([(df[!, flux_col_name] .- mean(df[!, flux_col_name])) ./ df[!, flux_err_col_name] for df in dfs]...) # Form the matrix whose columns are the whitened lightcurves
    sM = svd(M)

    (sM.U, sM.S) # Return the lightcurve basis and singular values
end

@doc raw"""
    full_detrend_basis_to_detrend_design_matrix(U, n_detrend)

Returns a design matrix suitable for detrending lightcurves composed of a
constant vector (to re-scale the mean lightcurve level) and `n_detrend` basis
vectors that have the largest signular values obtained from
`full_detrend_basis_and_singular_values`.  The constant vector is the first
column of the design matrix.
"""
function full_detrend_basis_to_detrend_design_matrix(U, n_detrend)
    nlc = size(U, 1)
    constant_vector = ones(nlc) / sqrt(nlc) # Normalized to unit norm

    [constant_vector U[:, 1:n_detrend]]
end

@doc raw"""
    detrend!(df, basis; flux_col_name="FLUX", flux_err_col_name="FLUX_ERR_EST")

Adducts columns `FLUX_DETREND`, `FLUX_DETREND_NORM`, and
`FLUX_DETREND_NORM_ERR_EST` to the data frame containing the detrended, and
normalized detrended flux and the normalized detrended flux uncertainty.
"""
function detrend!(df, basis; flux_col_name="FLUX", flux_err_col_name="FLUX_ERR_EST")
    x = basis \ df.FLUX # Least-squares solution.

    df[!, flux_col_name * "_DETREND"] = df[!, flux_col_name] - basis[:, 2:end] * x[2:end] # Subtract everything but the mean 
    df[!, flux_col_name * "_DETREND_NORM"] = df[!, flux_col_name * "_DETREND"] ./ (basis[:, 1] .* x[1]) # Normalize by the mean
    df[!, flux_col_name * "_DETREND_NORM_ERR_EST"] = df[!, flux_err_col_name] ./ (basis[:, 1] .* x[1]) # Normalize by the mean
end

@doc raw"""
    num_cbvs_chi2_threshold(U, dfs; threshold=2.0, flux_col_name="FLUX", flux_err_col_name="FLUX_ERR_EST")

Returns the number of CBVs to use for detrending based on the median delta-chi2
over stars.

Finds the first `n` where the median delta-chi2 from subtracting the `n+1`th CBV
is smaller than `threshold`.
"""
function num_cbvs_chi2_threshold(U, dfs; threshold=2.0, flux_col_name="FLUX", flux_err_col_name="FLUX_ERR_EST")
    for j in axes(U, 2)
        vhat = U[:,j]
        delta_chi2s = [(vhat' * df[!, flux_col_name] / df[1, flux_err_col_name]) .^ 2 for df in dfs]
        if median(delta_chi2s) < threshold
            return j - 1 # The first j where the median delta-chi2 is less than threshold
        end
    end
    return size(U, 2)
end

end # module KeplerDetrend

