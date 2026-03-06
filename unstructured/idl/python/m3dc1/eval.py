from __future__ import annotations

import numpy as np


def eval(field: np.ndarray, localpos: np.ndarray, theta: float, elm: int, operation: int = 1) -> float:
    """Port of eval.pro polynomial basis evaluation."""
    f = np.asarray(field)
    threed = f.shape[0] == 80

    mi = np.array([0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 3, 2, 1, 0], dtype=int)
    ni = np.array([0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 2, 3, 4, 5], dtype=int)

    mi1 = np.maximum(mi - 1, 0)
    mi2 = np.maximum(mi - 2, 0)
    mi3 = np.maximum(mi - 3, 0)
    ni1 = np.maximum(ni - 1, 0)
    ni2 = np.maximum(ni - 2, 0)
    ni3 = np.maximum(ni - 3, 0)

    lp = np.asarray(localpos, dtype=float)
    lp0 = np.array([1.0, lp[0], lp[0] ** 2, lp[0] ** 3, lp[0] ** 4, lp[0] ** 5], dtype=float)
    lp1 = np.array([1.0, lp[1], lp[1] ** 2, lp[1] ** 3, lp[1] ** 4, lp[1] ** 5], dtype=float)

    co = np.cos(theta)
    sn = np.sin(theta)

    op = int(operation)
    op2 = (op - 1) // 10
    op1 = op - op2 * 10

    if op1 == 1:
        temp = lp0[mi] * lp1[ni]
    elif op1 == 2:
        temp = co * mi * lp0[mi1] * lp1[ni] - sn * ni * lp0[mi] * lp1[ni1]
    elif op1 == 3:
        temp = sn * mi * lp0[mi1] * lp1[ni] + co * ni * lp0[mi] * lp1[ni1]
    elif op1 == 4:
        temp = (
            co * co * mi * mi1 * lp0[mi2] * lp1[ni]
            + sn * sn * ni * ni1 * lp0[mi] * lp1[ni2]
            - 2.0 * co * sn * ni * mi * lp0[mi1] * lp1[ni1]
        )
    elif op1 == 5:
        temp = (
            co * sn * mi * mi1 * lp0[mi2] * lp1[ni]
            - co * sn * ni * ni1 * lp0[mi] * lp1[ni2]
            + (co * co - sn * sn) * ni * mi * lp0[mi1] * lp1[ni1]
        )
    elif op1 == 6:
        temp = (
            sn * sn * mi * mi1 * lp0[mi2] * lp1[ni]
            + co * co * ni * ni1 * lp0[mi] * lp1[ni2]
            + 2.0 * co * sn * ni * mi * lp0[mi1] * lp1[ni1]
        )
    elif op1 == 7:
        temp = mi * mi1 * lp0[mi2] * lp1[ni] + ni * ni1 * lp1[ni2] * lp0[mi]
    elif op1 == 8:
        temp = (
            co * mi * mi1 * mi2 * lp0[mi3] * lp1[ni]
            - sn * ni * ni1 * ni2 * lp1[ni3] * lp0[mi]
            - sn * mi * mi1 * lp0[mi2] * ni * lp1[ni1]
            + co * mi * lp0[mi1] * ni * ni1 * lp1[ni2]
        )
    elif op1 == 9:
        temp = (
            sn * mi * mi1 * mi2 * lp0[mi3] * lp1[ni]
            + co * ni * ni1 * ni2 * lp1[ni3] * lp0[mi]
            + co * mi * mi1 * lp0[mi2] * ni * lp1[ni1]
            + sn * mi * lp0[mi1] * ni * ni1 * lp1[ni2]
        )
    else:
        temp = lp0[mi] * lp1[ni]

    if op2 == 0:
        summ = f[0:20, elm] * temp
        if threed:
            summ = summ + temp * (
                f[20:40, elm] * lp[2]
                + f[40:60, elm] * lp[2] ** 2
                + f[60:80, elm] * lp[2] ** 3
            )
    elif op2 == 1:
        if threed:
            summ = temp * (f[20:40, elm] + f[40:60, elm] * lp[2] * 2.0 + f[60:80, elm] * lp[2] ** 2 * 3.0)
        else:
            summ = 0.0
    elif op2 == 2:
        if threed:
            summ = temp * (f[40:60, elm] * 2.0 + f[60:80, elm] * lp[2] * 6.0)
        else:
            summ = 0.0
    else:
        summ = 0.0

    return float(np.sum(summ))
