"""Octree-based eigenvector reconstruction with interpolation.

Port of ``testOctreeEigReconstInterp.m``:

1. Build a saliency map: sum of normalised gradient magnitudes of all
   attribute volumes, smoothed with a ball kernel.
2. Zero-pad to a power-of-two cube and build an octree on the saliency map.
3. Run local PCA at the centre of every leaf node, or at its 8 corners when
   the node is large and its saliency varies.
4. Reconstruct per-voxel principal-component volumes, interpolating the 8
   corner eigenvectors where available, and write them as NRRD files.

All voxel positions are 1-based (as in MATLAB); array indexing subtracts 1.
"""
import os
import time

import numpy as np
from scipy import ndimage
from scipy.signal import fftconvolve

from eig_interp import create_bc_weights, eig_interp_3d
from local_pca import local_pca_space_dom_batch
from nrrd_io import nrrd_writer
from vol_octree import VolOctree
from volume_utils import aabb_intersect, load_volumes, show_montage

DEFAULT_FILE_LIST = ['./Toothuc0.20.nrrd', './Toothucderiv10.20.nrrd', './Toothucderiv20.20.nrrd']

# corner offsets (in units of the node size) for the 8-value nodes, in the
# order expected by eig_interp_3d
CORNER_OFFSETS = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                           [0, 1, 1], [1, 1, 1], [1, 0, 1], [0, 0, 1]])


def gradient_magnitude(vol):
    """Sobel gradient magnitude (``imgradient3`` default) with replicate borders."""
    g = [ndimage.sobel(vol, axis=a, mode='nearest') for a in range(3)]
    return np.sqrt(g[0] ** 2 + g[1] ** 2 + g[2] ** 2)


def ball_kernel(radius):
    """Normalised solid-ball averaging kernel of size (2*radius+1)^3."""
    r = np.arange(-radius, radius + 1)
    x, y, z = np.meshgrid(r, r, r, indexing='ij')
    k = (np.sqrt(x ** 2 + y ** 2 + z ** 2) <= radius).astype(np.float64)
    cnt = k.sum()
    return np.zeros_like(k) if cnt == 0 else k / cnt


def pad_offsets(shape, next_pow2):
    """Leading padding per axis to centre ``shape`` inside a ``next_pow2`` cube."""
    first = np.zeros(3, dtype=np.int64)
    for i in range(3):
        if shape[i] % 2 == 0:
            first[i] = (next_pow2 - shape[i]) // 2
        else:
            first[i] = next_pow2 // 2 - shape[i] // 2
    return first


def pad_volume(vol, next_pow2, pad_first):
    vp = np.zeros((next_pow2,) * 3)
    sx, sy, sz = vol.shape
    px, py, pz = pad_first
    vp[px:px + sx, py:py + sy, pz:pz + sz] = vol
    return vp


def test_octree_eig_reconst_interp(file_list=None, val_diff=0.01, do_sec_pc=True,
                                   octree_val_diff=0.1, R=3, show_plots=False,
                                   output_dir='.', seed=None):
    """Compute octree-reconstructed principal-component volumes.

    Returns the list of major-eigenvector volumes (one per attribute).

    ``val_diff`` is kept for signature compatibility; as in the MATLAB code the
    octree actually uses ``octree_val_diff`` (hard-coded 0.1 there).
    """
    if file_list is None:
        file_list = DEFAULT_FILE_LIST
    rng = np.random.default_rng(seed)
    if show_plots:
        import matplotlib.pyplot as plt

    # ---- load data volumes ------------------------------------------------
    vlist, n_vols, _ = load_volumes(file_list)
    vol_size = vlist[0].shape

    # gradient magnitude volume combining all attributes
    vgrad_sum = np.zeros(vol_size)
    for v in vlist:
        gmag = gradient_magnitude(v)
        gmin, gmax = gmag.min(), gmag.max()
        gmag = (gmag - gmin) / (gmax - gmin)
        if show_plots:
            show_montage(gmag)
        vgrad_sum += gmag

    # convolve with a ball whose radius is the diameter of the local neighbourhood
    rk = 2 * R + 1
    kernel = ball_kernel(rk)
    # kernel is symmetric, so convolution == MATLAB imfilter correlation (zero padding)
    vgrad_sum = fftconvolve(vgrad_sum, kernel, mode='same')
    if show_plots:
        plt.figure()
        plt.hist(vgrad_sum.ravel(), bins=100)
        show_montage(vgrad_sum, 'Filtered Gradient magnitude')

    # write out the saliency map
    nrrd_writer(os.path.join(output_dir, 'saliencyMap.nrrd'), vgrad_sum.astype(np.float32),
                [1, 1, 1], [0, 0, 0], 'raw')

    # ---- octree ------------------------------------------------------------
    t0 = time.time()
    print('Building the octree.')
    odim = np.array(vol_size)
    next_pow2 = int(2 ** np.ceil(np.log2(max(vol_size))))
    pad_first = pad_offsets(vol_size, next_pow2)
    vp = pad_volume(vgrad_sum, next_pow2, pad_first)
    if show_plots:
        show_montage(vp)
    ot = VolOctree(vp, [next_pow2] * 3, bin_val_diff=octree_val_diff)
    print('Elapsed time is %.6f seconds.' % (time.time() - t0))
    if show_plots:
        ot.plot(1, 1, 6)

    t0 = time.time()
    print('creating padded volumes')
    vpad_list = [pad_volume(v, next_pow2, pad_first) for v in vlist]
    del vlist
    print('Elapsed time is %.6f seconds.' % (time.time() - t0))

    # ---- local PCA on octree leaves ----------------------------------------
    org_min = 1 + pad_first
    org_max = odim + pad_first
    in_node_val_thres = 0.01
    n_vals_in_node = 8

    print('Calculating local PCA for octree nodes.')
    print('# of all nodes')
    print(ot.bin_count)
    t0 = time.time()

    leaf_ids = []        # octree node id per processed leaf
    leaf_n_vals = []     # 1 or 8
    pca_positions = []   # positions where local PCA is evaluated
    for i in range(ot.bin_count):
        if not ot.bin_is_leaf[i]:
            continue
        bounds = ot.bin_boundaries[i]
        bound_size = bounds[3:] - bounds[:3] + 1
        vmin, vmax = ot.bin_val_min_max[i]
        eight = bound_size[0] > 2 and vmax - vmin >= in_node_val_thres
        if not aabb_intersect(bounds[:3], bounds[3:], org_min, org_max):
            continue
        if eight:
            pca_positions.extend(bounds[:3] + CORNER_OFFSETS * bound_size)
            leaf_n_vals.append(n_vals_in_node)
        else:
            pca_positions.append(np.floor(0.5 * (bounds[3:] + bounds[:3])))
            leaf_n_vals.append(1)
        leaf_ids.append(i)

    eig_vecs, _ = local_pca_space_dom_batch(np.array(pca_positions), vpad_list, R, rng=rng)

    n_leaves = len(leaf_ids)
    maj_pc = np.zeros((n_leaves, n_vals_in_node, n_vols))
    maj_n = np.array(leaf_n_vals)
    sec_pc = np.zeros((n_leaves, n_vals_in_node, n_vols))
    sec_n = np.zeros(n_leaves, dtype=np.int64)
    k = 0
    for c in range(n_leaves):
        if maj_n[c] == n_vals_in_node:
            maj_pc[c] = eig_vecs[k:k + n_vals_in_node, :, 0]
            # NOTE: mirrors the MATLAB code: every corner gets the secondary PC of
            # the last corner, and the secondary count is left at 0.
            sec_pc[c] = eig_vecs[k + n_vals_in_node - 1, :, 1]
            k += n_vals_in_node
        else:
            maj_pc[c, 0] = eig_vecs[k, :, 0]
            sec_pc[c, 0] = eig_vecs[k, :, 1]
            sec_n[c] = 1
            k += 1

    print('# of leaf nodes')
    print(n_leaves)
    print('# of multivalue nodes')
    print(int(np.sum(maj_n > 1)))
    print('Elapsed time is %.6f seconds.' % (time.time() - t0))

    # ---- reconstruct PCA volumes from the octree ----------------------------
    pca_vols = np.zeros((n_vols,) + tuple(odim))
    pca_pad_vols = np.zeros((n_vols,) + (next_pow2,) * 3)
    sec_pca_vols = np.zeros((n_vols,) + tuple(odim)) if do_sec_pc else None

    min_level, max_level = 1, np.inf
    e2tri = np.zeros((n_vals_in_node, n_vols))  # persists across nodes, as in MATLAB
    for c, node_id in enumerate(leaf_ids):
        if ot.bin_depths[node_id] < min_level or ot.bin_depths[node_id] > max_level:
            continue
        bounds = ot.bin_boundaries[node_id]
        bound_size = bounds[3:] - bounds[:3] + 1
        # part of the node that lies inside the original (unpadded) volume
        lo = np.maximum(bounds[:3], org_min)
        hi = np.minimum(bounds[3:], org_max)
        if np.any(lo > hi):
            continue
        ps = tuple(slice(lo[a] - 1, hi[a]) for a in range(3))                    # padded
        os_ = tuple(slice(lo[a] - 1 - pad_first[a], hi[a] - pad_first[a]) for a in range(3))
        box = tuple(hi - lo + 1)

        if maj_n[c] > 1:
            etri = maj_pc[c]
            if do_sec_pc and sec_n[c] > 1:
                e2tri = sec_pc[c]
            axes = [np.arange(lo[a], hi[a] + 1) for a in range(3)]
            grid = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, 3)
            bc = create_bc_weights((grid - bounds[:3]) / bound_size)
            ev = eig_interp_3d(etri, bc).T.reshape((n_vols,) + box)
            pca_pad_vols[(slice(None),) + ps] = ev
            pca_vols[(slice(None),) + os_] = ev
            if do_sec_pc:
                sec_pca_vols[(slice(None),) + os_] = eig_interp_3d(e2tri, bc).T.reshape(
                    (n_vols,) + box)
            continue

        pca_info = maj_pc[c, 0]
        pca_pad_vols[(slice(None),) + ps] = pca_info[:, None, None, None]
        pca_vols[(slice(None),) + os_] = pca_info[:, None, None, None]
        if do_sec_pc:
            sec_pca_vols[(slice(None),) + os_] = sec_pc[c, 0][:, None, None, None]

    # ---- write out PCA volumes ---------------------------------------------
    for i in range(n_vols):
        name = os.path.splitext(os.path.basename(file_list[i]))[0]
        nrrd_writer(os.path.join(output_dir, '%s_Oct_majEigVolReconst.nrrd' % name),
                    pca_vols[i].astype(np.float32), [1, 1, 1], [0, 0, 0], 'raw')
        nrrd_writer(os.path.join(output_dir, '%s_PadOct_majEigVolReconst.nrrd' % name),
                    pca_pad_vols[i].astype(np.float32), [1, 1, 1], [0, 0, 0], 'raw')
        if do_sec_pc:
            nrrd_writer(os.path.join(output_dir, '%s_Oct_secEigVolReconst.nrrd' % name),
                        sec_pca_vols[i].astype(np.float32), [1, 1, 1], [0, 0, 0], 'raw')

    return [pca_vols[i] for i in range(n_vols)]
