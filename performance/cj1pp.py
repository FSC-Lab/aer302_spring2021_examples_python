#!/usr/bin/env python3

# cj1pp.m - power required and available
#
# created on: 25-Sep-00
# updated on: 27-Sep-00
#

## basic data for CP1, P45 lecture 08
from math import pi, sqrt

import matplotlib.pyplot as plt
import numpy as np

## basic data for CP1, P16 lecture 08
import data
from utils import stdatm


def main():
    ## power required at sea-level
    V0 = 152.4  # velocity for minimum power required
    C_L = data.W / (0.5 * data.rho_s * V0 ** 2 * data.S)
    C_D = data.C_D_0 + data.K * C_L ** 2
    D0 = (0.5 * data.rho_s * V0 ** 2 * data.S) * C_D

    # because V0 is the velocity for minimum power required, RR0=PR1
    PR0 = D0 * V0
    PR1 = data.W * sqrt(2 * data.WS / data.rho_s) / (C_L ** (3 / 2) / C_D)

    ## power required, power available plt.plots
    # at sea level

    Vs = np.arange(10, 410)
    qsS = 0.5 * data.rho_s * Vs ** 2 * data.S
    C_Ls = data.W / qsS
    C_Ds = data.C_D_0 + data.K * C_Ls ** 2
    Ds = qsS * C_Ds
    Fs = data.F_s
    PR_s = Ds * Vs
    PA_s = Fs * Vs

    # altitude h=6.7km
    h = 6705.6
    T_h, p_h, rho_h = stdatm(h)

    Vh = Vs * sqrt(data.rho_s / rho_h)
    qS = 0.5 * rho_h * Vh ** 2 * data.S
    C_L = data.W / qS
    C_D = data.C_D_0 + data.K * C_L ** 2
    Dh = qS * C_D
    Fh = (rho_h / data.rho_s) * Fs
    PR_h = Dh * Vh
    PA_h = PA_s * sqrt(rho_h / data.rho_s)

    # minimum power
    PRmin_s, iMin_s = (min(PR_s), np.argmin(PR_s))
    VPRmin_s = Vs[iMin_s]
    PRmin_h, iMin_h = (min(PR_h), np.argmin(PR_h))
    VPRmin_h = Vh[iMin_h]

    # minimum thrust
    TRmin_s, iMin_s = (min(Ds), np.argmin(Ds))
    VTRmin_s = Vs[iMin_s]

    # power required and power available at different altitude
    plt.figure(1)
    plt.plot(Vs, PR_s, "-", linewidth=2, label="PR at sea level")
    plt.plot(Vs, PA_s, "-", linewidth=2, label="PA at sea level")
    plt.plot(
        Vh, PR_h, "--", linewidth=2, label="PR at altitude {} m".format(h)
    )
    plt.plot(
        Vh, PA_h, "--", linewidth=2, label="PA at altitude {} m".format(h)
    )
    plt.plot(VPRmin_s, PRmin_s, "o", VPRmin_h, PRmin_h, "*")
    plt.grid()
    plt.legend()
    plt.title("Power Required and Available")
    plt.xlabel("velocity (m/s)")
    plt.ylabel("power (N m/s)")

    # power required and thrust available at sea-level
    plt.figure(2)
    plt.plot(
        Vs,
        PR_s / 100,
        "-",
        linewidth=2,
        label="PR at sea level (scaled by 100)",
    )
    plt.plot(Vs, Ds, "--", linewidth=2, label="TR at sea level")
    plt.plot(VPRmin_s, PRmin_s / 100, "o", VTRmin_s, TRmin_s, "*")
    plt.grid()
    plt.legend()
    plt.title("Power Required and Thrust Required at sea leval")
    plt.xlabel("velocity (m/s)")
    plt.ylabel("power (N m/s) or thrust (N)")

    plt.show()


if __name__ == "__main__":
    main()