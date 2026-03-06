from __future__ import annotations

from pathlib import Path

import numpy as np


def _read_ascii(path: Path) -> np.ndarray:
    return np.atleast_2d(np.loadtxt(str(path), dtype=float))


def read_coil_data(*, directory=".", rmp=False):
    d = Path(directory)
    coil_file = d / ("rmp_coil.dat" if rmp else "coil.dat")
    curr_file = d / ("rmp_current.dat" if rmp else "current.dat")
    if not coil_file.exists() or not curr_file.exists():
        return None
    coil = _read_ascii(coil_file)
    curr = _read_ascii(curr_file)
    n = min(coil.shape[0], curr.shape[0])
    if n == 0:
        return None
    out = np.zeros((10, n), dtype=float)
    out[0, :] = curr[:n, 0]
    out[1, :] = curr[:n, 1] if curr.shape[1] > 1 else curr[:n, 0]
    ccols = min(6, max(coil.shape[1] - 3, 0))
    if ccols > 0:
        out[2 : 2 + ccols, :] = coil[:n, 3 : 3 + ccols].T
    if ccols <= 6:
        out[8, :] = 1.0
        out[9, :] = 1.0
    return out

