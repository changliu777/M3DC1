from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

from .convert_units import convert_units
from .dimensions import dimensions
from .get_colors import get_colors
from .get_normalizations import get_normalizations
from .hdf5_file_test import hdf5_file_test
from .plot_legend import plot_legend
from .read_scalar import read_scalar


def _smooth_1d(x: np.ndarray, n: int) -> np.ndarray:
    if n <= 1:
        return np.asarray(x, dtype=float)
    k = np.ones(int(n), dtype=float) / float(n)
    return np.convolve(np.asarray(x, dtype=float), k, mode="same")


def _find_group_case_insensitive(h5: h5py.File, name: str):
    n = name.lower()
    for k in h5.keys():
        if k.lower() == n:
            return h5[k]
    return None


def _first_dataset(obj) -> np.ndarray | None:
    if isinstance(obj, h5py.Dataset):
        return np.asarray(obj[()])
    if isinstance(obj, h5py.Group):
        for k in ("_data", "data", "value"):
            if k in obj and isinstance(obj[k], h5py.Dataset):
                return np.asarray(obj[k][()])
        for k in ("KEHARMONICS", "BHARMONICS"):
            if k in obj:
                out = _first_dataset(obj[k])
                if out is not None:
                    return out
        for k in obj.keys():
            out = _first_dataset(obj[k])
            if out is not None:
                return out
    return None


def _read_harmonics(filename: str | Path, magnetic: bool) -> np.ndarray:
    key = "bharmonics" if magnetic else "keharmonics"
    with h5py.File(str(filename), "r") as h5:
        obj = _find_group_case_insensitive(h5, key)
        if obj is None:
            raise KeyError(f"Missing '{key}' in file.")
        arr = _first_dataset(obj)
        if arr is None:
            raise KeyError(f"Could not parse dataset for '{key}'.")
    out = np.asarray(arr, dtype=float)
    if out.ndim != 2:
        raise ValueError(f"Expected 2D harmonics array for '{key}', got {out.shape}.")
    return out


def _orient_harmonics(h: np.ndarray, nt: int) -> np.ndarray:
    arr = np.asarray(h, dtype=float)
    if arr.shape[1] == nt:
        return arr
    if arr.shape[0] == nt:
        return arr.T
    d0 = abs(arr.shape[0] - nt)
    d1 = abs(arr.shape[1] - nt)
    return arr if d1 <= d0 else arr.T


def plot_hmn(
    *,
    filename: str | Path = "C1.h5",
    maxn: int | None = None,
    growth: bool = False,
    outfile: str | Path | None = None,
    yrange=None,
    smooth: int | None = None,
    overplot: bool = False,
    thick: float | None = None,
    linestyle="-",
    ke: bool = False,
    me: bool = False,
    xscale: float = 1.0,
    labelx: float = 0.5,
    nolegend: bool = False,
    **kwargs,
):
    """Python port of plot_hmn.pro."""
    del kwargs, ke
    if not hdf5_file_test(filename):
        return None, None

    magnetic = bool(me)
    name = "Magnetic Energy" if magnetic else "Kinetic Energy"
    try:
        kehmn = _read_harmonics(filename, magnetic=magnetic)
    except Exception as exc:
        print(f"Error reading harmonics: {exc}")
        return None, None

    tmeta = read_scalar("time", filename=filename, return_meta=True)
    time = np.asarray(tmeta.data, dtype=float).reshape(-1)
    finite = np.isfinite(time)
    time = time[finite]
    if time.size == 0:
        return None, None
    kehmn = _orient_harmonics(kehmn, nt=time.size)
    ntimes = min(int(kehmn.shape[1]), int(time.size))
    kehmn = kehmn[:, :ntimes]
    time = time[:ntimes]

    d = dimensions(energy=1)
    b0, n0, l0, mi = get_normalizations(filename=filename)
    kehmn = convert_units(kehmn, d, b0=b0, n0=n0, l0=l0, mi=mi, filename=filename)

    if maxn is None:
        maxn = int(kehmn.shape[0])
    maxn = int(max(1, min(maxn, int(kehmn.shape[0]))))

    if outfile is not None:
        out = np.column_stack([time.reshape(-1), kehmn[:maxn, :].T])
        np.savetxt(str(outfile), out, fmt="%16.6e")

    keplot = np.asarray(kehmn[:maxn, :], dtype=float)
    tiny = np.finfo(float).tiny
    grate = np.zeros_like(keplot, dtype=float)
    for n in range(maxn):
        grate[n, :] = np.gradient(np.log(np.maximum(np.abs(keplot[n, :]), tiny)), time)

    tmp = grate if growth else keplot
    ytitle = "Growth Rate" if growth else name
    xtitle = f"t ({tmeta.units})" if tmeta.units else "t"

    if smooth is not None:
        for n in range(maxn):
            tmp[n, :] = _smooth_1d(tmp[n, :], int(smooth))

    if yrange is None:
        yrange = [float(np.nanmin(tmp)), float(np.nanmax(tmp))]

    colors = list(get_colors(maxn + 2))[1 : maxn + 1]
    if not overplot:
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.set_xlabel(xtitle)
        ax.set_ylabel(ytitle)
        if np.isfinite(yrange[0]) and np.isfinite(yrange[1]) and yrange[0] != yrange[1]:
            ax.set_ylim(float(yrange[0]), float(yrange[1]))
    else:
        ax = plt.gca()
        fig = ax.figure

    lw = float(thick) if thick is not None else None
    tscaled = np.asarray(time, dtype=float) * float(xscale)
    for n in range(maxn):
        c = colors[n % len(colors)] if colors else None
        ax.plot(tscaled, tmp[n, :], linestyle=linestyle, linewidth=lw, color=c)
        m = min(int(ntimes * float(labelx)), ntimes - 1)
        ax.text(tscaled[m], tmp[n, m], f"{n}", color=c if c is not None else "black")

    if not nolegend:
        names = [f"n={i}" for i in range(maxn)]
        plot_legend(names, colors=colors, linestyles=[linestyle] * maxn)

    return fig, ax
