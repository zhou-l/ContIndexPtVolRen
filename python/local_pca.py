"""Spatial-domain local PCA (port of ``localPCAspaceDom.m``).

The neighbourhood of a position is Monte-Carlo sampled with ``max_cnt``
random offsets in [-R, R]^3; each sample is trilinearly interpolated from all
attribute volumes and PCA is run on the resulting ``max_cnt x nVols`` matrix.
As in the MATLAB code, samples falling outside the volume are kept as rows of
zeros.

PCA follows MATLAB ``pca``: centre the columns, SVD, and flip each component
so that its largest-magnitude coefficient is positive.
"""
import numpy as np

from volume_utils import vol_interp


def _sign_convention(coeff):
    """Make the largest |element| of each column positive (batched)."""
    idx = np.argmax(np.abs(coeff), axis=-2)                       # (..., k)
    vals = np.take_along_axis(coeff, idx[..., None, :], axis=-2)  # (..., 1, k)
    return coeff * np.sign(vals)


def _pca_coeff(data):
    """Principal component coefficients of ``data`` (..., n, p), columns sorted by variance."""
    xc = data - data.mean(axis=-2, keepdims=True)
    _, _, vt = np.linalg.svd(xc, full_matrices=False)
    return _sign_convention(np.swapaxes(vt, -1, -2))


def local_pca_space_dom_batch(positions, vols, R, max_cnt=100, rng=None, chunk=2048):
    """Local PCA at many 1-based positions at once.

    Returns ``(eig_vecs, eig_vals)`` with shapes ``(N, nVols, nVols)`` (each
    column a principal component) and ``(N, nVols)``.
    """
    n_vols = len(vols)
    positions = np.atleast_2d(np.asarray(positions, dtype=np.float64))
    n = positions.shape[0]
    if n_vols < 1 or n == 0:
        return np.zeros((n, 0, 0)), np.zeros((n, 0))
    rng = np.random.default_rng() if rng is None else rng
    dim = np.array(vols[0].shape)

    eig_vecs = np.empty((n, n_vols, n_vols))
    for s in range(0, n, chunk):
        pos = positions[s:s + chunk]
        m = pos.shape[0]
        offsets = 2.0 * R * rng.random((m, max_cnt, 3)) - R
        npos = pos[:, None, :] + offsets
        valid = np.all((npos >= 1) & (npos <= dim), axis=-1)       # (m, max_cnt)
        d = np.zeros((m, max_cnt, n_vols))
        pts = npos[valid]
        for k, v in enumerate(vols):
            d[valid, k] = vol_interp(v, pts)
        eig_vecs[s:s + m] = _pca_coeff(d)[:, :, :n_vols]
    eig_vals = np.linalg.norm(eig_vecs, axis=1)  # vecnorm of each column (as in MATLAB)
    return eig_vecs, eig_vals


def local_pca_space_dom(pos, vols, R, max_cnt=100, rng=None):
    """Local PCA at a single 1-based position. Returns ``(eig_vecs, eig_vals)``."""
    if len(vols) < 1:
        return np.zeros((0, 0)), np.zeros(0)
    vecs, vals = local_pca_space_dom_batch(np.asarray(pos)[None, :], vols, R, max_cnt, rng)
    return vecs[0], vals[0]
