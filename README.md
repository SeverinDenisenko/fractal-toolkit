# Toolkit for fractal time series analysis

Initialy created for my specialist degree thesis

## Usage

Instal `gfortran`, `fortran_stdlib` and `gnuplot` (optianaly). Build by running:

```
make
```

Run tests:

```
make test
```

## API Reference

### Module `precision`
Kind parameters: `sp` (single), `dp` (double), `wp` (= double).

### Module `constants`
Parameters `pi`, `twopi`.

### Module `checks`
| Procedure | Description |
|---|---|
| `check(condition [, msg, code, ierr])` | Returns `.true.` on failure; aborts if no `ierr` passed, otherwise stores code and message |
| `last_message` | Variable holding the last error message |

### Module `version`
| Procedure | Description |
|---|---|
| `ver(v)` | Writes library version string into `v` |

### Module `math`
| Function | Description |
|---|---|
| `log2(a)` | Base-2 logarithm of `a` |
| `logn(a, n)` | Base-`n` logarithm of `a` |

### Module `integers`
| Function | Description |
|---|---|
| `up2power(n [, ierr])` | Smallest power of 2 >= `n` |
| `is2power(n)` | `.true.` if `n` is a power of 2 |

### Module `complex`
| Function | Description |
|---|---|
| `zroots(n, nn)` | First `nn` complex `n`-th roots of unity |

### Module `fourier`
| Subroutine | Description |
|---|---|
| `fft1d(data [, ierr])` | In-place forward FFT of complex vector (length must be power of 2) |
| `ifft1d(data [, ierr])` | In-place inverse FFT of complex vector |
| `fft2d(data [, ierr])` / `ifft2d(data [, ierr])` | Forward/inverse FFT along rows of a 2D array |
| `rfft1d(data, zdata [, ierr])` | Real-to-complex FFT of real vector into half-spectrum `zdata` |
| `irfft1d(data, zdata [, ierr])` | Complex-to-real inverse FFT (inverse of `rfft1d`) |
| `sfft1d`, `sfft2d` | Signed variants (`isign = ±1`) used internally |

### Module `conv`
| Subroutine | Description |
|---|---|
| `conv1d_full(y, x, b [, ierr])` | Full linear convolution of `x` with kernel `b` via FFT (`size(y) = nx + nb - 1`) |
| `conv1d_same(y, x, b [, ierr])` | Convolution truncated to length of `x` |

### Module `io`
| Subroutine | Description |
|---|---|
| `write_table_file(name, array1 [, array2] [, ierr])` | Writes one or two columns of numbers to text file `name` |
| `read_single_column(name, array [, ierr])` | Reads single-column file into allocated array |

### Module `stat`
| Procedure | Description |
|---|---|
| `mean(S)` | Arithmetic mean |
| `variance(S [, i])` | Variance; optional correction `i` added to sample size (`-1` unbiased default, `0` biased) |
| `percentile(S, q, pct [, ierr])` | Linear-interpolated percentiles `q` (%) of data |
| `fd_bins(S)` | Optimal histogram bin count via Freedman-Diaconis rule |
| `skewness(S)` | Asymmetry of distribution (positive = longer right tail) |
| `kurtosis(S)` | Excess kurtosis vs normal distribution |
| `normal_cdf(x, mu, sigma)` | Normal cumulative distribution function (elemental) |
| `ks_statistic(S)` | Kolmogorov-Smirnov statistic vs fitted normal (input sorted ascending) |
| `cvm_criterion(S)` | Cramer-von Mises criterion vs fitted normal (input sorted ascending) |
| `fit_normal(S, mu, sigma [, se_mu, se_sigma, ks, cvm, sk, ku])` | Fits normal distribution and returns standard errors, goodness-of-fit statistics and moments |

### Module `autoreg`
| Procedure | Description |
|---|---|
| `yw_ar_coeff(phi, S [, ierr])` | AR coefficients via least-squares Yule-Walker equations |
| `complex_yw_ar_coeff(phi, S [, ierr])` | Yule-Walker AR coefficients for complex series `S` |
| `burg_ar_coeff(phi, S)` | AR coefficients via Burg method |
| `complex_burg_ar_coeff(phi, S)` | Burg AR coefficients for complex series `S` |
| `ar_predict(phi, S, i)` | One-step-ahead AR prediction of `S(i)` from previous values |
| `complex_ar_predict(phi, S, i)` | One-step-ahead AR prediction for complex data |
| `ar_freq_response(phi, f, H [, ierr])` | Frequency response `H(ω)` of the AR filter at normalized frequencies `f` |
| `complex_ar_freq_response(phi, f, H [, ierr])` | Frequency response of a complex-coefficient AR filter |

### Module `spectra`
| Procedure | Description |
|---|---|
| `psd_size(n)` | Number of PSD points for series of length `n` (`n/2 + 1`) |
| `freq_full(f, dt)` | Frequencies from 0 up to sampling frequency for time step `dt` |
| `freq_nyquist(f, dt)` | Frequencies from 0 up to Nyquist frequency for time step `dt` |
| `ar_psd(f, P, S, dt, phi [, ierr])` | Parametric power spectral density from given AR coefficients |
| `complex_ar_psd(f, P, S, dt, phi [, ierr])` | PSD from complex AR coefficients |
| `berg_psd(f, P, S, dt, m [, ierr])` | PSD estimated by Burg method of order `m` |
| `complex_berg_psd(f, P, S, dt, m [, ierr])` | Burg PSD for complex series `S` of order `m` |
| `yw_psd(f, P, S, dt, m [, ierr])` | PSD estimated by Yule-Walker method of order `m` |
| `complex_yw_psd(f, P, S, dt, m [, ierr])` | Yule-Walker PSD for complex series `S` of order `m` |

### Module `frac`
| Procedure | Description |
|---|---|
| `frac_diff_coeff(coeff, order)` | Maclaurin-series coefficients for fractional differencing of given order |
| `frac_diff(coeff, series_in, series_out [, ierr])` | Fractional difference using precomputed coefficients |
| `frac_diff_simple(order, series_in, series_out [, ierr])` | Fractional difference computing coefficients internally |

### Module `solvers`
| Subroutine | Description |
|---|---|
| `linregress(x, y, k, a [, sigma2, k_err, a_err, ierr])` | Ordinary least squares fit `y = k*x + a` with residuals variance and parameter errors |
| `powerregress(x, y, a, c [, sigma2, a_err, c_err, ierr])` | Power-law fit `y = c*x^a` (log-log least squares); requires positive data |

### Module `interp`
| Subroutine | Description |
|---|---|
| `cubic_resample(X, Y, X1, Y1 [, ierr])` | Cubic spline interpolation of `Y` evaluated at new points `X1` into `Y1` |
| `cubic_solve(X, Y, b, c, d [, ierr])` | Computes cubic spline coefficients `b`, `c`, `d` |
| `cubic_interp(X, Y, b, c, d, X1, Y1 [, ierr])` | Evaluates cubic spline at points `X1` using precomputed coefficients |

### Module `generators`
| Subroutine | Description |
|---|---|
| `generate_white(array [, mu, sigma, seed])` | Uniform white noise (reproducible via seed) |
| `generate_gauss(array [, mu, sigma, seed])` | Gaussian white noise (reproducible via seed) |
| `generate_gauss_integrate(series, intorder [, seed] [, ierr])` | Fractionally integrated Gaussian noise |
| `generate_white_integrate(series, intorder [, seed] [, ierr])` | Fractionally integrated uniform noise |
| `generate_color(series, a [, seed] [, ierr])` | Colored noise with PSD ~ 1/f^a |
| `generate_fgn(series, H [, seed] [, ierr])` | Exact fractional Gaussian noise with Hurst exponent `H` (Davies-Harte method) |

### Module `hurst`
| Procedure | Description |
|---|---|
| `slope_to_hurst(a)` | Hurst exponent from PSD slope `a` where PSD ~ 1/f^a |
| `estimate_hurst_psd(f, P, H, a, H_err, a_err, sigma2 [, ierr])` | Hurst estimate from precomputed power spectrum |
| `complex_estimate_hurst_psd(f, P, H, a, H_err, a_err, sigma2 [, ierr])` | Hurst estimate from precomputed complex spectrum |
| `estimate_hurst_berg(series, m, H, a, H_err, a_err, sigma2 [, ierr])` | Hurst exponent via Burg PSD of order `m` |
| `complex_estimate_hurst_berg(series, m, H, a, H_err, a_err, sigma2 [, ierr])` | Hurst exponent via Burg PSD for complex series |
| `estimate_hurst_yw(series, m, H, a, H_err, a_err, sigma2 [, ierr])` | Hurst exponent via Yule-Walker PSD of order `m` |
| `complex_estimate_hurst_yw(series, m, H, a, H_err, a_err, sigma2 [, ierr])` | Hurst exponent via Yule-Walker PSD for complex series |
| `rs_chart_size(n)` | Length of R/S chart arrays for series of size `n` |
| `rs_chart(series, RS, N [, ierr])` | R/S statistic chart over chunk sizes N = 2, 4, ... |
| `complex_rs_chart(series, RS, N [, ierr])` | R/S statistic chart for complex series (range = diameter of cumulative walk) |
| `estimate_hurst_rs(series, H, H_err, sigma2 [, ierr])` | Hurst exponent via rescaled range (R/S) analysis |
| `complex_estimate_hurst_rs(series, H, H_err, sigma2 [, ierr])` | Hurst exponent via R/S analysis for complex series |
| `estimate_hurst_lssd(X, p, q, H, H_err [, eps, maxiter] [, ierr])` | Hurst exponent via Koutsoyiannis Least Squares Standard Deviation estimator |

All procedures follow the same error convention: pass an optional `integer ierr` argument to receive an error code instead of program termination.

The non-prefixed procedures accept real series and are implemented as thin bindings that cast the input to complex (zero imaginary part) and call the corresponding `complex_` variant. Results for real input are identical to the underlying complex algorithm.
