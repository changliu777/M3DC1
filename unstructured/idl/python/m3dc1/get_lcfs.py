from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from .read_field import read_field
from .read_lcfs import read_lcfs


def get_lcfs(psi=None, x=None, z=None, *, psival=None, axis=None, filename="C1.h5", points=200, slice=0, **kwargs):
    """
    Return an approximate LCFS contour path as shape (2, N).
    """
    if psi is None or x is None or z is None:
        pmeta = read_field("psi", filename=filename, slices=slice, points=points, equilibrium=True, return_meta=True, **kwargs)
        psi2d = np.asarray(pmeta.data)[0, :, :] if np.asarray(pmeta.data).ndim == 3 else np.asarray(pmeta.data)
        xv = np.asarray(pmeta.r, dtype=float).reshape(-1)
        zv = np.asarray(pmeta.z, dtype=float).reshape(-1)
    else:
        arr = np.asarray(psi)
        psi2d = arr[0, :, :] if arr.ndim == 3 else arr
        xv = np.asarray(x, dtype=float).reshape(-1)
        zv = np.asarray(z, dtype=float).reshape(-1)

    if psival is None:
        lc = read_lcfs(filename=filename, slice=slice, return_meta=True, **kwargs)
        psival = float(lc.psilim)

    cs = plt.contour(xv, zv, psi2d.T, levels=[float(psival)])
    try:
        if not cs.collections or not cs.collections[0].get_paths():
            return np.zeros((2, 0), dtype=float)
        path = cs.collections[0].get_paths()[0].vertices
        return np.asarray(path.T, dtype=float)
    finally:
        plt.clf()

