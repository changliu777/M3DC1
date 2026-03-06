from __future__ import annotations

import inspect
import numpy as np

from .dimensions import dimensions
from .flux_average_field import flux_average_field
from .flux_coordinates import flux_coordinates
from .parse_units import parse_units
from .read_field import read_field

_READ_FIELD_KW = set(inspect.signature(read_field).parameters.keys())


def flux_average(
    field,
    *,
    psi=None,
    i0=None,
    x=None,
    z=None,
    t=None,
    r0=None,
    flux=None,
    nflux=None,
    area=None,
    bins=None,
    points=None,
    name=None,
    symbol=None,
    units=None,
    integrate=False,
    complex=False,
    abs=False,
    phase=False,
    stotal=False,
    fac=None,
    fc=None,
    elongation=None,
    filename="C1.h5",
    linear=False,
    linfac=1.0,
    mks=False,
    cgs=False,
    return_meta=False,
    **kwargs,
):
    """
    Python port of flux_average.pro (core behavior used by plotting routines).
    """
    read_kwargs = {k: v for k, v in kwargs.items() if k in _READ_FIELD_KW}

    if fc is None:
        if points is None:
            points = 200
        nb = int(bins if bins is not None else points)
        if psi is None or x is None or z is None:
            p = read_field(
                "psi",
                filename=filename,
                slices=int(t or 0),
                points=int(points),
                equilibrium=True,
                return_meta=True,
                cgs=cgs,
                mks=mks,
                **read_kwargs,
            )
            psi = np.asarray(p.data)[0, :, :] if np.asarray(p.data).ndim == 3 else np.asarray(p.data)
            x = np.asarray(p.r, dtype=float).reshape(-1)
            z = np.asarray(p.z, dtype=float).reshape(-1)
        fc = flux_coordinates(psi0=psi, i0=i0, x=x, z=z, points=int(points), fbins=nb, tbins=int(points), filename=filename, **kwargs)
    if flux is not None:
        flux[:] = np.asarray(fc.psi)
    if nflux is not None:
        nflux[:] = np.asarray(fc.psi_norm)
    if area is not None:
        area[:] = np.asarray(fc.area)

    if isinstance(field, str):
        fkey = field.strip().lower()
        if fkey in {"safety factor", "q"}:
            vals = np.abs(np.asarray(fc.q, dtype=float))
            title = "Safety Factor"
            sym = "$q$"
            u = ""
            if return_meta:
                return vals, title, sym, u, fc
            return vals
        if fkey == "rho":
            vals = np.asarray(fc.rho, dtype=float)
            title = "rho"
            sym = "$\\rho$"
            u = parse_units(dimensions(l0=1), cgs=cgs, mks=mks)
            if return_meta:
                return vals, title, sym, u, fc
            return vals
        if fkey == "flux_t":
            vals = np.asarray(fc.flux_tor, dtype=float)
            if return_meta:
                return vals, "Toroidal Flux", "$\\Phi_t$", "", fc
            return vals
        if fkey == "flux_p":
            vals = np.asarray(fc.flux_pol, dtype=float)
            if return_meta:
                return vals, "Poloidal Flux", "$\\Phi_p$", "", fc
            return vals
        if fkey == "volume":
            vals = np.asarray(fc.V, dtype=float)
            u = parse_units(dimensions(l0=3), cgs=cgs, mks=mks)
            if return_meta:
                return vals, "Volume", "$V$", u, fc
            return vals

        meta = read_field(
            field,
            slices=int(t or 0),
            points=int(points or fc.m),
            linear=linear,
            complex=complex,
            abs=abs,
            phase=phase,
            fac=fac,
            linfac=linfac,
            filename=filename,
            cgs=cgs,
            mks=mks,
            return_meta=True,
            **read_kwargs,
        )
        arr = np.asarray(meta.data)
        farr = arr[0, :, :] if arr.ndim == 3 else arr
        vals = flux_average_field(
            farr,
            psi=psi,
            x=np.asarray(x),
            z=np.asarray(z),
            fc=fc,
            bins=bins,
            integrate=integrate,
            surface_weight=stotal,
            filename=filename,
            **kwargs,
        )
        title = str(meta.symbol)
        sym = str(meta.symbol) if not integrate and not stotal else str(meta.symbol)
        u = str(meta.units)
        if return_meta:
            return np.asarray(vals), title, sym, u, fc
        return np.asarray(vals)

    vals = flux_average_field(
        field,
        psi=psi,
        x=np.asarray(x),
        z=np.asarray(z),
        fc=fc,
        bins=bins,
        integrate=integrate,
        surface_weight=stotal,
        filename=filename,
        **kwargs,
    )
    if return_meta:
        return np.asarray(vals), "", "", "", fc
    return np.asarray(vals)
