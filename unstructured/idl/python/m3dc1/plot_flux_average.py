from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np

from .dimensions import dimensions
from .flux_at_q import flux_at_q
from .flux_average import flux_average
from .get_colors import get_colors
from .get_slice_time import get_slice_time
from .parse_units import parse_units
from .plot_legend import plot_legend
from .read_parameter import read_parameter


def _smooth_1d(x: np.ndarray, n: int) -> np.ndarray:
    if n <= 1:
        return np.asarray(x)
    k = np.ones(int(n), dtype=float) / float(n)
    return np.convolve(np.asarray(x, dtype=float), k, mode="same")


def _write_outfile(outfile: str | Path, flux: np.ndarray, fa: np.ndarray, complex_mode: bool = False) -> None:
    f = np.asarray(flux, dtype=float).reshape(-1)
    y = np.asarray(fa)
    if complex_mode:
        arr = np.column_stack([f, np.real(y).reshape(-1), np.imag(y).reshape(-1)])
    else:
        arr = np.column_stack([f, np.real(y).reshape(-1)])
    np.savetxt(str(outfile), arr, fmt="%16.6e")


def plot_flux_average(
    field,
    time=0,
    *,
    filename: str | Path | Sequence[str | Path] = "C1.h5",
    complex: bool = False,
    color=None,
    names=None,
    bins: int | None = None,
    linear: bool = False,
    xlog: bool = False,
    ylog: bool = False,
    overplot: bool = False,
    fac: float = 1.0,
    lcfs: bool = False,
    normalized_flux: bool = True,
    points: int = 200,
    linestyle="-",
    linfac: float = 1.0,
    sum: bool = False,
    minor_radius: bool = False,
    smooth: int | None = None,
    t=None,
    rms: bool = False,
    bw: bool = False,
    srnorm: bool = False,
    last: bool = False,
    mks: bool = False,
    cgs: bool = False,
    q_contours=None,
    rho: bool = False,
    integrate: bool = False,
    multiply_flux: float = 0.0,
    abs: bool = False,
    phase: bool = False,
    stotal: bool = False,
    nolegend: bool = False,
    outfile: str | Path | None = None,
    val_at_q=None,
    flux_at_q_out=None,
    regularize: bool = False,
    **kwargs,
):
    """
    Python port of plot_flux_average.pro.
    """
    if last:
        time = int(read_parameter("ntime", filename=filename)) - 1

    if isinstance(field, (list, tuple, np.ndarray)) and not isinstance(field, str):
        ff = list(field)
        n = len(ff)
        cols = [color] * n if color is not None else (["black"] * n if bw else list(get_colors(n + 2))[1 : n + 1])
        lss = list(linestyle) if isinstance(linestyle, (list, tuple)) else [linestyle] * n
        lf = list(linfac) if isinstance(linfac, (list, tuple, np.ndarray)) else [linfac] * n
        for i, f in enumerate(ff):
            qcon = q_contours if i == 0 else None
            plot_flux_average(
                f,
                time,
                filename=filename,
                complex=complex,
                color=cols[i],
                bins=bins,
                linear=linear,
                xlog=xlog,
                ylog=ylog,
                overplot=(i > 0) or overplot,
                fac=fac,
                lcfs=lcfs,
                normalized_flux=normalized_flux,
                points=points,
                linestyle=lss[i],
                linfac=lf[i],
                sum=sum,
                minor_radius=minor_radius,
                smooth=smooth,
                t=t,
                rms=rms,
                bw=bw,
                srnorm=srnorm,
                last=last,
                mks=mks,
                cgs=cgs,
                q_contours=qcon,
                rho=rho,
                integrate=integrate,
                multiply_flux=multiply_flux,
                abs=abs,
                phase=phase,
                stotal=stotal,
                nolegend=True,
                regularize=regularize,
                **kwargs,
            )
        if names is not None and not nolegend:
            plot_legend(list(names), colors=cols, linestyles=lss)
        if plt.get_fignums():
            ax = plt.gca()
            return ax.figure, ax
        return None, None

    if isinstance(filename, (list, tuple)) and (len(filename) > 1) and not sum:
        fl = list(filename)
        n = len(fl)
        cols = [color] * n if color is not None else (["black"] * n if bw else list(get_colors(n + 2))[1 : n + 1])
        lss = list(linestyle) if isinstance(linestyle, (list, tuple)) else [linestyle] * n
        tlist = list(time) if isinstance(time, (list, tuple, np.ndarray)) else [time] * n
        for i, fn in enumerate(fl):
            plot_flux_average(
                field,
                tlist[i],
                filename=fn,
                complex=complex,
                color=cols[i],
                bins=bins,
                linear=linear,
                xlog=xlog,
                ylog=ylog,
                overplot=(i > 0) or overplot,
                fac=fac,
                lcfs=lcfs,
                normalized_flux=normalized_flux,
                points=points,
                linestyle=lss[i],
                linfac=linfac,
                sum=sum,
                minor_radius=minor_radius,
                smooth=smooth,
                t=t,
                rms=rms,
                bw=bw,
                srnorm=srnorm,
                last=last,
                mks=mks,
                cgs=cgs,
                q_contours=q_contours,
                rho=rho,
                integrate=integrate,
                multiply_flux=multiply_flux,
                abs=abs,
                phase=phase,
                stotal=stotal,
                nolegend=True,
                outfile=None if not isinstance(outfile, (list, tuple)) else outfile[i],
                regularize=regularize,
                **kwargs,
            )
        if names is not None and not nolegend:
            plot_legend(list(names), colors=cols, linestyles=lss)
        if plt.get_fignums():
            ax = plt.gca()
            return ax.figure, ax
        return None, None

    if isinstance(time, (list, tuple, np.ndarray)) and np.asarray(time).size > 1:
        tv = np.asarray(time).reshape(-1)
        n = tv.size
        cols = [color] * n if color is not None else (["black"] * n if bw else list(get_colors(n + 2))[1 : n + 1])
        lss = list(linestyle) if isinstance(linestyle, (list, tuple)) else [linestyle] * n
        labels = []
        for i in range(n):
            plot_flux_average(
                field,
                int(tv[i]),
                filename=filename,
                complex=complex,
                color=cols[i],
                bins=bins,
                linear=linear,
                xlog=xlog,
                ylog=ylog,
                overplot=(i > 0) or overplot,
                fac=fac,
                lcfs=lcfs,
                normalized_flux=normalized_flux,
                points=points,
                linestyle=lss[i],
                linfac=linfac,
                sum=sum,
                minor_radius=minor_radius,
                smooth=smooth,
                t=t,
                rms=rms,
                bw=bw,
                srnorm=srnorm,
                last=last,
                mks=mks,
                cgs=cgs,
                q_contours=q_contours,
                rho=rho,
                integrate=integrate,
                multiply_flux=multiply_flux,
                abs=abs,
                phase=phase,
                stotal=stotal,
                nolegend=True,
                regularize=regularize,
                **kwargs,
            )
            ts = float(np.asarray(get_slice_time(filename=filename, slice=int(tv[i]), cgs=cgs, mks=mks)).reshape(-1)[0])
            tu = parse_units(dimensions(t0=1), cgs=cgs, mks=mks)
            labels.append(f"t = {ts:g} {tu}")
        if names is None:
            names = labels
        if names is not None and not nolegend:
            plot_legend(list(names), colors=cols, linestyles=lss)
        if plt.get_fignums():
            ax = plt.gca()
            return ax.figure, ax
        return None, None

    vals, title, symbol, units, fc = flux_average(
        field,
        t=time,
        bins=bins,
        points=points,
        linear=linear,
        fac=fac,
        linfac=linfac,
        filename=filename,
        integrate=integrate,
        complex=complex,
        abs=abs,
        phase=phase,
        stotal=stotal,
        mks=mks,
        cgs=cgs,
        return_meta=True,
        **kwargs,
    )
    fa = np.asarray(vals)
    if fa.size <= 1:
        return None, None

    ytitle = symbol
    if units:
        ytitle = f"{ytitle} ({units})"
    xtitle = "psi"

    flux = np.asarray(fc.psi, dtype=float)
    nflux = np.asarray(fc.psi_norm, dtype=float)
    lcfs_psi = float(nflux[-1]) if nflux.size > 0 else 1.0

    if rms:
        fa2 = np.asarray(
            flux_average(
                np.asarray(read_field(field, filename=filename, slices=int(time), points=points, linear=linear, complex=complex, abs=abs, phase=phase, return_meta=False, **kwargs))
                ** 2,
                t=time,
                bins=bins,
                points=points,
                filename=filename,
                return_meta=False,
                **kwargs,
            )
        )
        fa = np.sqrt(np.maximum(1.0 - np.asarray(fa) ** 2 / np.maximum(fa2, np.finfo(float).tiny), 0.0))
        title = f"Poloidal Deviation of {title}"

    if multiply_flux:
        flux = flux * float(multiply_flux)
        nflux = nflux * float(multiply_flux)

    if srnorm:
        flux = np.sqrt(np.maximum(nflux, 0.0))
        xtitle = "sqrt(psi_n)"
        lcfs_psi = 1.0
    elif minor_radius:
        rvals = flux_average("r", t=time, bins=bins, points=points, linear=linear, fac=fac, filename=filename, mks=mks, cgs=cgs, return_meta=False, **kwargs)
        flux = np.asarray(rvals, dtype=float)
        xtitle = "minor radius"
    elif rho:
        flux = np.asarray(fc.rho, dtype=float)
        xtitle = "rho"
        lcfs_psi = 1.0
    elif normalized_flux:
        flux = nflux
        xtitle = "psi_n"
        lcfs_psi = 1.0

    if regularize:
        denom = np.sum(fa)
        if abs(denom) > np.finfo(float).tiny:
            fa = fa / denom

    if smooth is not None:
        fa = _smooth_1d(np.asarray(fa, dtype=float), int(smooth))

    if not overplot:
        plt.figure(figsize=(7, 4))
    plt.plot(flux, np.real(fa), color=color, linestyle=linestyle)
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    if title:
        plt.title(title)
    if xlog:
        plt.xscale("log")
    if ylog:
        plt.yscale("log")

    if lcfs:
        y0, y1 = plt.ylim()
        plt.plot([lcfs_psi, lcfs_psi], [y0, y1], linestyle="--", color=color if color is not None else "k")

    if q_contours is not None:
        fvals = flux_at_q(q_contours, points=points, filename=filename, normalized_flux=bool(normalized_flux), **kwargs)
        fv = np.asarray(fvals, dtype=float).reshape(-1)
        if fv.size > 0 and fv[0] > 0:
            y0, y1 = plt.ylim()
            for xq in fv:
                plt.plot([xq, xq], [y0, y1], linestyle="-", color=color if color is not None else "k", linewidth=0.7)
            valsq = np.interp(fv, np.asarray(flux, dtype=float), np.asarray(np.real(fa), dtype=float))
            if val_at_q is not None:
                val_at_q = valsq
            if flux_at_q_out is not None:
                flux_at_q_out = fv

    if outfile is not None:
        _write_outfile(outfile, flux, fa, complex_mode=complex)

    ax = plt.gca()
    return ax.figure, ax
