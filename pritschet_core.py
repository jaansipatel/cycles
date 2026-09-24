"""Signal-processing functions for the cycles pipeline, matched to the code
behind Pritschet and Santander 2020 (tsantander/PritschetSantander2020_NI_Hormones).

Four places where their code differs from what the paper says, kept on
purpose because the code is what produced the published results:
  1. the frequency band is isolated with wavelets, not a bandpass filter
  2. intensity scaling divides by the median, not the mean
  3. the detrending step hits the noise regressors too, not just the data
  4. the motion regressors are [R, R(t-1), R^2, R(t-1)^2], no intercept
"""

import numpy as np
import pywt


def modwt_filters(wavelet="db8"):
    # pywt stores the filters back to front, so flip them, then rescale by
    # 1/sqrt(2), which is what makes this the MODWT rather than the plain DWT
    w = pywt.Wavelet(wavelet)
    g = np.asarray(w.dec_lo)[::-1] / np.sqrt(2.0)
    h = np.asarray(w.dec_hi)[::-1] / np.sqrt(2.0)
    return g, h


def modwt(x, n_levels, wavelet="db8"):
    """Wavelet transform that treats the signal as a loop. Returns
    (details [n_levels x N], smooth [N]).

    The wrap-around uses np.roll. The old version glued one extra copy of
    the signal on as padding, which fails once the stretched filter gets
    longer than the signal itself (level 7 on short runs); roll wraps as
    many times as needed, so any length and level combination is safe.
    """
    x = np.asarray(x, dtype=float)
    n = x.shape[0]
    g, h = modwt_filters(wavelet)
    v = x.copy()
    details = np.empty((n_levels, n))
    for j in range(1, n_levels + 1):
        s = 2 ** (j - 1)
        wj = np.zeros(n)
        vj = np.zeros(n)
        for l, (gl, hl) in enumerate(zip(g, h)):
            rolled = np.roll(v, s * l)
            wj += hl * rolled
            vj += gl * rolled
        details[j - 1] = wj
        v = vj
    return details, v


def imodwt(details, smooth, wavelet="db8"):
    # exact inverse of modwt() above
    g, h = modwt_filters(wavelet)
    v = smooth.copy()
    n_levels = details.shape[0]
    for j in range(n_levels, 0, -1):
        s = 2 ** (j - 1)
        v_new = np.zeros_like(v)
        for l, (gl, hl) in enumerate(zip(g, h)):
            v_new += hl * np.roll(details[j - 1], -s * l) + gl * np.roll(v, -s * l)
        v = v_new
    return v


def modwt_band(x, keep_scales, wavelet="db8"):
    """Rebuild the signal from the chosen scales only, dropping the rest.
    Works column by column on 2D input (time x channels)."""
    x = np.asarray(x, dtype=float)
    if x.ndim == 2:
        return np.column_stack([modwt_band(x[:, k], keep_scales, wavelet)
                                for k in range(x.shape[1])])
    n_levels = max(keep_scales)
    details, smooth = modwt(x, n_levels, wavelet)
    mask = np.zeros((n_levels, 1))
    mask[[j - 1 for j in keep_scales]] = 1.0
    return imodwt(details * mask, np.zeros_like(smooth), wavelet)


def scale_band_hz(scale, tr):
    # each scale covers one octave: scale j spans fs/2^(j+1) to fs/2^j
    fs = 1.0 / tr
    return fs / 2 ** (scale + 1), fs / 2 ** scale


def solve_scales(tr, f_lo, f_hi, max_scale=10):
    """Which scales overlap the target band [f_lo, f_hi] Hz at this TR."""
    keep = [j for j in range(1, max_scale + 1)
            if scale_band_hz(j, tr)[0] < f_hi and scale_band_hz(j, tr)[1] > f_lo]
    if not keep:
        raise ValueError(f"no MODWT scale overlaps [{f_lo}, {f_hi}] Hz at TR={tr}")
    return keep


def grand_median_scale(data, mask, target=1000.0):
    """Scale a run so the median value inside the mask hits the target.
    Returns (scaled data, factor); reuse the factor on any copy of the same
    run that has to stay on the same scale."""
    med = np.median(data[mask])
    factor = target / med
    return data * factor, factor


def friston24_ar1(rp):
    """Expand the 6 motion estimates into 24 regressors: the values, the
    same values one step back, and both squared. No intercept column; the
    detrending step already removes means."""
    rp = np.asarray(rp, dtype=float)
    lag = np.vstack([np.zeros((1, rp.shape[1])), rp[:-1]])
    return np.hstack([rp, lag, rp ** 2, lag ** 2])


def detrend_projector(n):
    # matrix that removes the mean and any straight-line drift
    p = np.column_stack([np.ones(n), np.arange(n, dtype=float)])
    return np.eye(n) - p @ np.linalg.pinv(p)


def nuisance_regress(y, x):
    """Remove the noise regressors x from y. Both sides get detrended first
    (difference 3 in the list above); returns what is left of y."""
    m = detrend_projector(y.shape[0])
    my, mx = m @ y, m @ x
    beta, *_ = np.linalg.lstsq(mx, my, rcond=None)
    return my - mx @ beta


def first_eigenvariate(y):
    """One summary signal for a region (time x voxels), the SPM way: the
    strongest shared component, flipped if needed so it tracks the region
    mean rather than its mirror image."""
    yc = y - y.mean(axis=0)
    u, s, vt = np.linalg.svd(yc, full_matrices=False)
    e = u[:, 0] * s[0] / np.sqrt(y.shape[1])
    if np.corrcoef(e, y.mean(axis=1))[0, 1] < 0:
        e = -e
    return e
