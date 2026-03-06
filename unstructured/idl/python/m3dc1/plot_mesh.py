from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from .make_label import make_label
from .read_parameter import read_parameter
from .read_mesh import read_mesh


def plot_mesh(*, mesh=None, oplot=False, boundary=False, iso: bool = False, filename="C1.h5", slice=0, **kwargs):
    """
    Plot triangular mesh edges.
    """
    m = read_mesh(filename=filename, slice=slice) if mesh is None else mesh
    el = np.asarray(m.elements, dtype=float)
    if oplot:
        ax = plt.gca()
        fig = ax.figure
    else:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_xlabel(make_label("R", l0=1, **kwargs))
        ax.set_ylabel(make_label("Z", l0=1, **kwargs))

    if el.ndim != 2 or el.shape[0] < 6:
        if iso:
            ax.set_aspect("equal", adjustable="box")
        return fig, ax

    try:
        version = read_parameter("version", filename=filename)
        print(f"Output version = {version}")
    except Exception:
        pass

    nelms = int(m.nelms)
    minr = np.array([np.inf, np.inf], dtype=float)
    maxr = np.array([-np.inf, -np.inf], dtype=float)
    for i in range(nelms):
        a = float(el[0, i])
        b = float(el[1, i])
        c = float(el[2, i])
        t = float(el[3, i])
        x0 = float(el[4, i])
        y0 = float(el[5, i])
        bound = int(el[6, i]) if el.shape[0] > 6 else 0

        p1 = np.asarray([x0, y0], dtype=float)
        p2 = p1 + np.asarray([(b + a) * np.cos(t), (b + a) * np.sin(t)], dtype=float)
        p3 = p1 + np.asarray([b * np.cos(t) - c * np.sin(t), b * np.sin(t) + c * np.cos(t)], dtype=float)
        minr[0] = min(minr[0], p1[0], p2[0], p3[0])
        minr[1] = min(minr[1], p1[1], p2[1], p3[1])
        maxr[0] = max(maxr[0], p1[0], p2[0], p3[0])
        maxr[1] = max(maxr[1], p1[1], p2[1], p3[1])

        if not boundary:
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="0.35", linewidth=0.4)
            ax.plot([p2[0], p3[0]], [p2[1], p3[1]], color="0.35", linewidth=0.4)
            ax.plot([p3[0], p1[0]], [p3[1], p1[1]], color="0.35", linewidth=0.4)
            continue

        if bound & 1:
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="k", linewidth=1.1)
        if bound & 2:
            ax.plot([p2[0], p3[0]], [p2[1], p3[1]], color="k", linewidth=1.1)
        if bound & 4:
            ax.plot([p3[0], p1[0]], [p3[1], p1[1]], color="k", linewidth=1.1)

    if iso:
        ax.set_aspect("equal", adjustable="box")

    print(f"Elements: {nelms}")
    print(f"sqrt(nodes) (estimated): {np.sqrt(nelms / 2.0)}")
    if np.isfinite(minr).all() and np.isfinite(maxr).all():
        print(f"Width: {maxr[0] - minr[0]}")
        print(f"Height: {maxr[1] - minr[1]}")

    return fig, ax
