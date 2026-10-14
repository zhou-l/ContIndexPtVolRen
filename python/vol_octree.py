"""Volume octree (port of the parts of ``VolOctree.m`` used by the pipeline).

A cubic, power-of-two volume is recursively split into 8 equal children
while the value range (max - min) inside a bin exceeds ``bin_val_diff``.

Bin boundaries use MATLAB's 1-based, inclusive voxel indices
``[xmin ymin zmin xmax ymax zmax]``. Bin ids are 0-based; ``bin_parents``
holds -1 for the root.
"""
import numpy as np


class VolOctree:
    def __init__(self, vol_data, vol_dim, bin_val_diff=0.0, max_depth=np.inf,
                 max_size=np.inf, min_size=1, bbox_min=None, bbox_max=None):
        next_pow2 = int(2 ** np.ceil(np.log2(max(vol_dim))))
        self.vol = vol_data
        self.vol_dim = tuple(vol_dim)
        self.properties = dict(
            binValDiff=bin_val_diff, maxDepth=max_depth, maxSize=max_size, minSize=min_size,
            bboxMin=list(bbox_min) if bbox_min is not None else [1, 1, 1],
            bboxMax=list(bbox_max) if bbox_max is not None else [next_pow2] * 3)

        self._bounds = [self.properties['bboxMin'] + self.properties['bboxMax']]
        self._depths = [0]
        self._parents = [-1]
        self._minmax = [(float(np.min(vol_data)), float(np.max(vol_data)))]

        self._divide([0])

        self.bin_boundaries = np.array(self._bounds, dtype=np.int64)
        self.bin_depths = np.array(self._depths, dtype=np.int64)
        self.bin_parents = np.array(self._parents, dtype=np.int64)
        self.bin_val_min_max = np.array(self._minmax, dtype=np.float64)
        self.depth_max = int(self.bin_depths.max())
        has_children = np.zeros(self.bin_count, dtype=bool)
        has_children[self.bin_parents[self.bin_parents >= 0]] = True
        self.bin_is_leaf = ~has_children

    @property
    def bin_count(self):
        return len(self._bounds)

    def _divide(self, starting_bins):
        p = self.properties
        for bin_no in starting_bins:
            if self._depths[bin_no] + 1 >= p['maxDepth']:
                continue
            b = self._bounds[bin_no]
            edge = [b[3] - b[0] + 1, b[4] - b[1] + 1, b[5] - b[2] + 1]
            min_edge, max_edge = min(edge), max(edge)
            if min_edge < p['minSize']:
                continue
            if min_edge == p['minSize']:
                continue
            old_count = self.bin_count
            vmin, vmax = self._minmax[bin_no]
            if vmax - vmin > p['binValDiff']:
                self._divide_bin(bin_no)
                self._divide(range(old_count, self.bin_count))
                continue
            if max_edge > p['maxSize']:
                self._divide_bin(bin_no)
                self._divide(range(old_count, self.bin_count))
                continue

    def _divide_bin(self, bin_no):
        b = self._bounds[bin_no]
        old_min = np.array(b[:3])
        old_max = np.array(b[3:])
        new_edge = (old_max - old_min + 1) // 2
        new_div = old_min + new_edge - 1
        # minMidMax = [oldMin newDiv newDiv+1 oldMax], 0-based column indices below
        mmm = np.concatenate([old_min, new_div, new_div + 1, old_max])
        rows = [[0, 1, 2, 3, 4, 5],
                [0, 1, 8, 3, 4, 11],
                [0, 7, 2, 3, 10, 5],
                [0, 7, 8, 3, 10, 11],
                [6, 1, 2, 9, 4, 5],
                [6, 1, 8, 9, 4, 11],
                [6, 7, 2, 9, 10, 5],
                [6, 7, 8, 9, 10, 11]]
        depth = self._depths[bin_no] + 1
        for r in rows:
            nb = [int(v) for v in mmm[r]]
            self._bounds.append(nb)
            self._depths.append(depth)
            self._parents.append(bin_no)
            block = self.vol[nb[0] - 1:nb[3], nb[1] - 1:nb[4], nb[2] - 1:nb[5]]
            self._minmax.append((float(block.min()), float(block.max())))

    def plot(self, step_size=1, min_level=1, max_level=np.inf, ax=None, **plot_kwargs):
        """Draw the bounding boxes of leaf bins whose depth is in [min_level, max_level]."""
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        if ax is None:
            ax = plt.figure().add_subplot(projection='3d')
        edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
                 (0, 4), (1, 5), (2, 6), (3, 7)]
        segments = []
        for i in range(0, self.bin_count, step_size):
            if not self.bin_is_leaf[i]:
                continue
            lvl = self.bin_depths[i]
            if lvl < min_level or lvl > max_level:
                continue
            x0, y0, z0, x1, y1, z1 = self.bin_boundaries[i]
            c = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
                 (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
            segments.extend([c[a], c[b]] for a, b in edges)
        ax.add_collection3d(Line3DCollection(segments, colors=plot_kwargs.get('color', 'b'),
                                             linewidths=plot_kwargs.get('linewidth', 0.5)))
        lo, hi = self.properties['bboxMin'], self.properties['bboxMax']
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
        return ax
