"""Volume helpers: loading, interpolation, resizing, AABB test, plotting.

Ports of ``loadVolumes.m``, ``volInterp.m``, ``myVolDownSample.m`` and
``aabbIntersect.m``.

Positions passed to these functions follow MATLAB's 1-based convention so
that all index arithmetic in the ported algorithms stays identical.
"""
import numpy as np

from nrrd_io import nhdr_nrrd_read

NAN_FILL_VALUE = -6.0


def load_volumes(file_list):
    """Load NRRD/NHDR volumes, replacing NaN with -6.

    Returns ``(vlist, n_vols, vol_info_list)``.
    """
    vlist = []
    info_list = []
    for filename in file_list:
        header = nhdr_nrrd_read(filename, True)
        info_list.append(header)
        v = np.array(header['data'], dtype=np.float64)
        v[np.isnan(v)] = NAN_FILL_VALUE
        vlist.append(v)
    return vlist, len(file_list), info_list


def vol_interp(vol, pos):
    """Trilinear interpolation of ``vol`` at 1-based position(s) ``pos``.

    ``pos`` may be a single (3,) position or an (N, 3) array. Positions must
    lie in [1, dim] per axis.
    """
    pos = np.asarray(pos, dtype=np.float64)
    single = pos.ndim == 1
    pos = np.atleast_2d(pos)
    p0 = np.floor(pos).astype(np.int64)
    a = pos - p0
    dims = np.array(vol.shape)
    i0 = p0 - 1                                   # to 0-based
    # floor(pos)+1 only exceeds the volume when pos == dim exactly; its weight is 0 then
    i1 = np.minimum(i0 + 1, dims - 1)
    x0, y0, z0 = i0[:, 0], i0[:, 1], i0[:, 2]
    x1, y1, z1 = i1[:, 0], i1[:, 1], i1[:, 2]
    ax, ay, az = a[:, 0], a[:, 1], a[:, 2]

    c00 = (1 - ax) * vol[x0, y0, z0] + ax * vol[x1, y0, z0]
    c10 = (1 - ax) * vol[x0, y1, z0] + ax * vol[x1, y1, z0]
    c01 = (1 - ax) * vol[x0, y0, z1] + ax * vol[x1, y0, z1]
    c11 = (1 - ax) * vol[x0, y1, z1] + ax * vol[x1, y1, z1]
    c0 = (1 - ay) * c00 + ay * c10
    c1 = (1 - ay) * c01 + ay * c11
    ev = (1 - az) * c0 + az * c1
    return ev[0] if single else ev


def _cubic(x):
    ax = np.abs(x)
    ax2 = ax * ax
    ax3 = ax2 * ax
    return ((1.5 * ax3 - 2.5 * ax2 + 1) * (ax <= 1)
            + (-0.5 * ax3 + 2.5 * ax2 - 4 * ax + 2) * ((ax > 1) & (ax <= 2)))


def _contributions(in_len, out_len, scale, antialiasing=True):
    """Weights/indices for one axis, following MATLAB's imresize algorithm."""
    kernel_width = 4.0
    if scale < 1 and antialiasing:
        h = lambda x: scale * _cubic(scale * x)
        kernel_width = kernel_width / scale
    else:
        h = _cubic
    x = np.arange(1, out_len + 1, dtype=np.float64)[:, None]
    u = x / scale + 0.5 * (1 - 1 / scale)
    left = np.floor(u - kernel_width / 2)
    p = int(np.ceil(kernel_width)) + 2
    ind = left + np.arange(p)[None, :]
    weights = h(u - ind)
    weights = weights / weights.sum(axis=1, keepdims=True)
    aux = np.concatenate([np.arange(1, in_len + 1), np.arange(in_len, 0, -1)])
    ind = aux[np.mod(ind.astype(np.int64) - 1, aux.size)]
    keep = np.any(weights != 0, axis=0)
    return weights[:, keep], ind[:, keep] - 1


def _resize_axis(vol, axis, out_len):
    in_len = vol.shape[axis]
    if in_len == out_len:
        return vol
    w, ind = _contributions(in_len, out_len, out_len / in_len)
    v = np.moveaxis(vol, axis, 0)
    out = np.einsum('ok,ok...->o...', w, v[ind])
    return np.moveaxis(out, 0, axis)


def imresize3(vol, out_shape):
    """Cubic resize with antialiasing, mimicking MATLAB ``imresize3``."""
    out = np.asarray(vol, dtype=np.float64)
    for ax in range(3):
        out = _resize_axis(out, ax, int(out_shape[ax]))
    return out


def my_vol_downsample(vlist, sample_rate_per_axis):
    """Resample every volume by ``sample_rate_per_axis`` (``myVolDownSample.m``)."""
    if len(vlist) < 1:
        return []
    dim = np.array(vlist[0].shape)
    ds_dim = np.floor(dim * sample_rate_per_axis).astype(int)
    return [imresize3(v, ds_dim) for v in vlist]


def aabb_intersect(a_min, a_max, b_min, b_max):
    """True if axis-aligned boxes A and B overlap (inclusive bounds)."""
    return all(a_min[k] <= b_max[k] and a_max[k] >= b_min[k] for k in range(3))


def show_montage(vol, title=None):
    """Show all z-slices of a 3-D volume in a grid (like MATLAB ``montage``)."""
    import matplotlib.pyplot as plt
    nz = vol.shape[2]
    ncols = int(np.ceil(np.sqrt(nz)))
    nrows = int(np.ceil(nz / ncols))
    h, w = vol.shape[0], vol.shape[1]
    canvas = np.full((nrows * h, ncols * w), np.nan)
    for z in range(nz):
        r, c = divmod(z, ncols)
        canvas[r * h:(r + 1) * h, c * w:(c + 1) * w] = vol[:, :, z]
    plt.figure()
    plt.imshow(canvas, cmap='gray', vmin=np.nanmin(vol), vmax=np.nanmax(vol))
    plt.axis('off')
    if title:
        plt.title(title)
