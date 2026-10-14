"""Eigenvector interpolation inside an octree cell (port of ``eigInterp3D.m``).

Corner order of ``etri_list`` rows (same as the MATLAB code)::

       3-E2-2
      /|   /|
     0-E1-1 |
     | 4-E3-5
     |/   |/
     7-E4-6

i.e. 0:(0,0,0) 1:(1,0,0) 2:(1,1,0) 3:(0,1,0) 4:(0,1,1) 5:(1,1,1) 6:(1,0,1) 7:(0,0,1).
"""
import numpy as np


def create_bc_weights(ds_pos, ids_pos=(0, 0, 0)):
    """Per-axis weights ``[1-x, x, 1-y, y, 1-z, z]``.

    ``ds_pos`` may be (3,) or (N, 3); the result is (6,) or (6, N).
    """
    d = np.asarray(ds_pos, dtype=np.float64) - np.asarray(ids_pos, dtype=np.float64)
    d = d.T
    return np.stack([1.0 - d[0], d[0], 1.0 - d[1], d[1], 1.0 - d[2], d[2]])


def max_pairwise_angle(etri_list):
    """Largest angle between any two (unit) row vectors of ``etri_list``."""
    e = np.asarray(etri_list, dtype=np.float64)
    dots = e @ e.T
    iu = np.triu_indices(e.shape[0], k=1)
    # clip guards against |dot| > 1 from rounding (MATLAB's acos would go complex)
    return float(np.max(np.arccos(np.clip(dots[iu], -1.0, 1.0))))


def eig_interp_3d(etri_list, bc):
    """Spherical trilinear interpolation of 8 corner eigenvectors.

    ``bc`` is the (6,) or (6, N) output of :func:`create_bc_weights`.
    Returns (nVols,) or (N, nVols) unit vectors.
    """
    e = np.asarray(etri_list, dtype=np.float64)
    bc = np.asarray(bc, dtype=np.float64)
    single = bc.ndim == 1
    bc = bc.reshape(6, -1)
    theta = max_pairwise_angle(e)
    if theta == 0:
        out = np.repeat(e[0][None, :], bc.shape[1], axis=0)
        return out[0] if single else out

    w = np.sin(bc * theta) / np.sin(theta)                  # (6, N)
    w1, w2, w3, w4, w5, w6 = (wi[:, None] for wi in w)
    # Be careful with the order!
    ex1 = w1 * e[0] + w2 * e[1]
    ex2 = w2 * e[2] + w1 * e[3]
    ex3 = w1 * e[4] + w2 * e[5]
    ex4 = w2 * e[6] + w1 * e[7]
    ey1 = w3 * ex1 + w4 * ex2
    ey2 = w4 * ex3 + w3 * ex4
    ez = w5 * ey1 + w6 * ey2
    with np.errstate(invalid='ignore', divide='ignore'):
        ev = ez / np.linalg.norm(ez, axis=1, keepdims=True)
    return ev[0] if single else ev
