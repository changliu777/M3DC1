from __future__ import annotations

import re
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np

from .compensate_renorm import compensate_renorm
from .get_colors import get_colors
from .make_label import make_label
from .plot_legend import plot_legend
from .power_spectrum import power_spectrum
from .read_parameter import read_parameter
from .read_scalar import read_scalar


def _plain(s: str) -> str:
    return re.sub(r"![A-Za-z0-9]+", "", str(s)).strip()


def _smooth_1d(x: np.ndarray, n: int) -> np.ndarray:
    if n <= 1:
        return x
    k = np.ones(n, dtype=float) / n
    return np.convolve(x, k, mode="same")


def _write_outfile(outfile: str | Path, tdata: np.ndarray, data: np.ndarray) -> None:
    d = np.asarray(data)
    if d.ndim == 1:
        arr = np.column_stack([tdata, d])
    else:
        arr = np.column_stack([tdata, d.T])
    np.savetxt(str(outfile), arr, fmt="%16.6e")


def plot_scalar(
    scalarname: str,
    x: float | Sequence[float] | None = None,
    filename: str | Path | Sequence[str | Path] = "C1.h5",
    names: Sequence[str] | None = None,
    overplot: bool = False,
    difference: bool = False,
    ylog: bool = False,
    xlog: bool = False,
    absolute_value: bool = False,
    power_spectrum_on: bool = False,
    per_length: bool = False,
    growth_rate: bool = False,
    bw: bool = False,
    nolegend: bool = False,
    cgs: bool = False,
    mks: bool = False,
    linestyle: str | Sequence[str] | None = None,
    color=None,
    outfile: str | Path | None = None,
    smooth: int | None = None,
    compensate_renorm_on: bool = False,
    integrate: bool = False,
    xscale: float = 1.0,
    ipellet: int = 0,
    factor: float = 1.0,
    versus: str | None = None,
    xabs: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Python port of plot_scalar.pro.
    """
    if isinstance(filename, (list, tuple)):
        flist = list(filename)
        nfiles = len(flist)
        if names is None:
            names = [str(f) for f in flist]
        if isinstance(x, (int, float)):
            xvals = [float(x)] * nfiles
        elif isinstance(x, (list, tuple, np.ndarray)):
            xvals = list(x)
        else:
            xvals = [None] * nfiles

        if bw:
            ls = list(linestyle) if isinstance(linestyle, (list, tuple)) else ["-"] * nfiles
            cs = list(color) if isinstance(color, (list, tuple)) else ["black"] * nfiles
        else:
            ls = list(linestyle) if isinstance(linestyle, (list, tuple)) else ["-"] * nfiles
            cs = list(color) if isinstance(color, (list, tuple)) else get_colors(nfiles)

        last = None
        for i, f in enumerate(flist):
            last = plot_scalar(
                scalarname,
                x=xvals[i],
                filename=f,
                overplot=(i > 0) or overplot,
                difference=difference,
                ylog=ylog,
                xlog=xlog,
                absolute_value=absolute_value,
                power_spectrum_on=power_spectrum_on,
                per_length=per_length,
                growth_rate=growth_rate,
                bw=bw,
                nolegend=True,
                cgs=cgs,
                mks=mks,
                linestyle=ls[i],
                color=cs[i],
                outfile=None,
                smooth=smooth,
                compensate_renorm_on=compensate_renorm_on,
                integrate=integrate,
                xscale=xscale,
                ipellet=ipellet,
                factor=factor,
                versus=versus,
                xabs=xabs,
            )
        if not nolegend and names:
            plot_legend(list(names), linestyles=ls, colors=cs, nolegend=nolegend)
        if last is None:
            return np.array([]), np.array([])
        return last

    meta = read_scalar(
        scalarname,
        filename=filename,
        integrate=integrate,
        growth_rate=growth_rate,
        ipellet=ipellet,
        cgs=cgs,
        mks=mks,
        return_meta=True,
    )
    data = np.asarray(meta.data, dtype=float) * factor
    time = np.asarray(meta.time, dtype=float)
    title = _plain(meta.title) if meta.title else scalarname
    ytitle = _plain(meta.symbol) if meta.symbol else scalarname
    units = meta.units

    if compensate_renorm_on:
        if data.ndim == 1:
            data = compensate_renorm(data)
        else:
            for i in range(data.shape[0]):
                data[i, :] = compensate_renorm(data[i, :])

    if data.size <= 1:
        return np.array([]), np.array([])

    if power_spectrum_on:
        if data.ndim == 1:
            data, tdata = power_spectrum(data, t=float(np.max(time)))
        else:
            out = np.zeros_like(data)
            tdata = None
            for i in range(data.shape[0]):
                out[i, :], freq = power_spectrum(data[i, :], t=float(np.max(time)))
                tdata = freq
            data = out
            assert tdata is not None
        xtitle = make_label("Frequency", t0=-1, cgs=cgs, mks=mks)
    else:
        tdata = time
        xtitle = make_label("Time", t0=1, cgs=cgs, mks=mks)
        print(xtitle)

    if versus:
        vmeta = read_scalar(versus, filename=filename, ipellet=ipellet, cgs=cgs, mks=mks, return_meta=True)
        tdata = np.asarray(vmeta.data, dtype=float)
        if xabs:
            tdata = np.abs(tdata)
        vlabel = _plain(vmeta.symbol) if vmeta.symbol else versus
        xtitle = f"{vlabel} ({vmeta.units})" if vmeta.units else vlabel

    if per_length:
        itor = int(read_parameter("itor", filename=filename))
        if itor == 1:
            rzero = float(read_parameter("rzero", filename=filename, cgs=cgs, mks=mks))
            data = data / rzero

    if difference:
        if data.ndim == 1:
            data = data - data[0]
        else:
            data = data - data[:, [0]]
        ytitle = "Delta " + ytitle

    if absolute_value:
        data = np.abs(data)
        ytitle = "|" + ytitle + "|"

    if growth_rate:
        ytitle = make_label("gamma", t0=-1, cgs=cgs, mks=mks)
    elif units:
        ytitle = f"{ytitle} ({units})"

    if smooth is not None and smooth > 1:
        if data.ndim == 1:
            data = _smooth_1d(data, smooth)
        else:
            for i in range(data.shape[0]):
                data[i, :] = _smooth_1d(data[i, :], smooth)

    if outfile is not None:
        _write_outfile(outfile, np.asarray(tdata), np.asarray(data))

    if x is None:
        if not overplot:
            plt.figure(figsize=(10, 5))
        if data.ndim == 1:
            plt.plot(tdata * xscale, data, linestyle=linestyle or "-", color=color)
        else:
            cs = get_colors(data.shape[0] + 3)
            for i in range(data.shape[0]):
                c = color if color is not None else cs[i + 1]
                plt.plot(tdata * xscale, data[i, :], linestyle=linestyle or "-", color=c)
    else:
        xv = float(x)
        z = float(np.asarray(data).reshape(-1)[-1])
        if not overplot:
            plt.figure(figsize=(6, 4))
        plt.plot([xv], [z], linestyle=linestyle or "None", marker="o", color=color)

    plt.title(title)
    plt.xlabel(_plain(xtitle))
    plt.ylabel(_plain(ytitle))
    if xlog:
        plt.xscale("log")
    if ylog:
        plt.yscale("log")
    plt.tight_layout()

    return np.asarray(tdata), np.asarray(data)
