"""Gold-standard per-voxel eigenvector computation.

Ports of ``myEigInfoCompute.m`` and ``testGoldEigReconstruct.m``: local PCA
is evaluated at every voxel of the (optionally resampled) volumes.
"""
import os
import time

import numpy as np

from local_pca import local_pca_space_dom_batch
from nrrd_io import nrrd_writer
from volume_utils import load_volumes, my_vol_downsample


def my_eig_info_compute(file_list, sample_rate_per_axis, R=3, rng=None):
    """Returns ``(maj_eig_vols, sec_eig_vols, maj_eig_val_vols, full_dim, ds_dim)``."""
    vlist, n_vols, _ = load_volumes(file_list)
    full_dim = np.array(vlist[0].shape)

    vs_list = my_vol_downsample(vlist, sample_rate_per_axis)
    ds_dim = np.array(vs_list[0].shape)

    # every voxel, in MATLAB's x-fastest order, 1-based
    xs, ys, zs = (np.arange(1, d + 1) for d in ds_dim)
    zz, yy, xx = np.meshgrid(zs, ys, xs, indexing='ij')
    positions = np.stack([xx.ravel(), yy.ravel(), zz.ravel()], axis=1)

    eig_vecs, eig_vals = local_pca_space_dom_batch(positions, vs_list, R, rng=rng)

    shape = tuple(ds_dim)
    maj_eig_vols = [np.zeros(shape) for _ in range(n_vols)]
    sec_eig_vols = [np.zeros(shape) for _ in range(n_vols)]
    idx = (positions[:, 0] - 1, positions[:, 1] - 1, positions[:, 2] - 1)
    for d in range(n_vols):
        maj_eig_vols[d][idx] = eig_vecs[:, d, 0]   # the major eigenvector
        sec_eig_vols[d][idx] = eig_vecs[:, d, 1]   # the second eigenvector
    maj_eig_val_vols = np.zeros(shape)
    maj_eig_val_vols[idx] = eig_vals[:, 0]
    return maj_eig_vols, sec_eig_vols, maj_eig_val_vols, full_dim, ds_dim


def test_gold_eig_reconstruct(file_list, sample_rate_per_axis=1, output_dir='.', seed=None):
    """Compute and write the gold-standard major/secondary eigenvector volumes."""
    print('Computing the gold standard.')
    t0 = time.time()
    rng = np.random.default_rng(seed)
    maj, sec, _, _, _ = my_eig_info_compute(file_list, sample_rate_per_axis, rng=rng)
    print('Elapsed time is %.6f seconds.' % (time.time() - t0))
    for i in range(len(maj)):
        name = os.path.splitext(os.path.basename(file_list[i]))[0]
        fn = '%s_eigMaj_gold%i_%0.2f.nrrd' % (name, i + 1, sample_rate_per_axis)
        nrrd_writer(os.path.join(output_dir, fn), maj[i].astype(np.float32),
                    [1, 1, 1], [0, 0, 0], 'raw')
        fn2 = '%s_eigSec_gold%i_%0.2f.nrrd' % (name, i + 1, sample_rate_per_axis)
        nrrd_writer(os.path.join(output_dir, fn2), sec[i].astype(np.float32),
                    [1, 1, 1], [0, 0, 0], 'raw')
    return maj
