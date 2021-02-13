#!/usr/bin/env python3

# cp1pp.m - power required and available, including climbing performance
#
# created on: 10.Feb.2020
# updated on:
#

from math import pi, sqrt

import matplotlib.pyplot as plt
import numpy as np

## basic data for CP1, P16 lecture 08
from data import cp1data
from utils import stdatm


def main():
    ## power required at sea-level P30-32 lecture 08
    V0 = 60.96  # velocity for minimum power required
    C_L = cp1data.W / (0.5 * cp1data.rho_s * V0 ** 2 * cp1data.S)
    C_D = cp1data.C_D_0 + cp1data.K * C_L ** 2
    D0 = (0.5 * cp1data.rho_s * V0 ** 2 * cp1data.S) * C_D

    # because V0 is the velocity for minimum power required, RR0=PR1
    PR0 = D0 * V0
    PR1 = (
        cp1data.W
        * sqrt(2 * cp1data.WS / cp1data.rho_s)
        / (C_L ** (3 / 2) / C_D)
    )

    ## power required and power available plt.plots
    # at sea level

    Vs = np.arange(10, 210)
    qsS = 0.5 * cp1data.rho_s * Vs ** 2 * cp1data.S
    C_Ls = cp1data.W / qsS
    C_Ds = cp1data.C_D_0 + cp1data.K * C_Ls ** 2
    Ds = qsS * C_Ds
    PR_s = Ds * Vs
    PA_s = cp1data.PAs
    RC_s = (PA_s - PR_s) / cp1data.W  # climb rate
    gamma_s = np.arcsin(RC_s / Vs) * 180 / pi  # flight path angle

    # at h=6km
    h = 6000
    [T_h, p_h, rho_h] = stdatm(h)
    PAh = cp1data.PAs * sqrt(rho_h / cp1data.rho_s)
    Vh = Vs * sqrt(cp1data.rho_s / rho_h)
    Vh = Vs
    qS = 0.5 * rho_h * Vh ** 2 * cp1data.S
    C_L = cp1data.W / qS
    C_D = cp1data.C_D_0 + cp1data.K * C_L ** 2
    Dh = qS * C_D
    PR_h = Dh * Vh
    PA_h = PAh
    RC_h = (PA_h - PR_h) / cp1data.W
    gamma_h = np.arcsin(RC_h / Vs) * 180 / pi

    # minimum power from graphics
    [PRmin_s, iMin_s] = (min(PR_s), np.argmin(PR_s))
    VPRmin_s = 10 + iMin_s * 1
    [PRmin_h, iMin_h] = (min(PR_h), np.argmin(PR_h))
    VPRmin_h = (10 + iMin_h * 1) * sqrt(cp1data.rho_s / rho_h)

    # minimum power from calculation (analytical)
    VPRmin_s_calc = sqrt(
        2 * cp1data.WS / cp1data.rho_s / sqrt(3 * cp1data.C_D_0 / cp1data.K)
    )
    VPRmin_h_calc = sqrt(
        2 * cp1data.WS / rho_h / sqrt(3 * cp1data.C_D_0 / cp1data.K)
    )

    # plt.plot
    plt.figure(1)
    plt.plot(Vs, PR_s, "-", label="PR at sea level", linewidth=2)
    plt.axhline(y=PA_s, linestyle="-", label="PA at sea level", linewidth=2)
    plt.plot(
        Vh, PR_h, "--", label="PR at altitude {} m".format(h), linewidth=2
    )
    plt.axhline(
        y=PA_h,
        linestyle="--",
        label="PA at altitude {} m".format(h),
        linewidth=2,
    )
    plt.plot(VPRmin_s, PRmin_s, "o", VPRmin_h, PRmin_h, "*")
    plt.grid()
    plt.axis([0, 250, 0, 2e06])
    plt.title("Power Required and Available for CP1")
    plt.xlabel("velocity (m/cp1data.s)")
    plt.ylabel("power (N m/cp1data.s)")

    ## hodograph (approximation)
    Vv_s = RC_s
    Vh_s = Vs * np.cos(gamma_s * pi / 180)
    Vv_h = RC_h
    Vh_h = Vh * np.cos(gamma_h * pi / 180)

    plt.figure(3)
    plt.plot(Vh_s, Vv_s, "-", label="sea-level", linewidth=2)
    plt.plot(Vh_h, Vv_h, "--", label="6000m", linewidth=2)
    plt.grid()
    plt.axis([0, 100, 0, 10])
    plt.title("Hodograph for [CP1]")
    plt.xlabel("Vh (m/cp1data.s)")
    plt.ylabel("Vv (m/cp1data.s)")

    ## other performance computations
    # maximum rate of climb
    VRC_max_s_calc = VPRmin_s_calc  # sea level
    VRC_max_h_calc = VPRmin_h_calc  # h=6km
    C_D_rcmax = 4 * cp1data.C_D_0  # cp1data.C_D_0 = 1/3 cp1data.K CL**2
    RCmax_s = (cp1data.PAs - PRmin_s) / cp1data.W  # sea level
    RCmax_h = (
        cp1data.PAs * sqrt(rho_h / cp1data.rho_s) - PRmin_h
    ) / cp1data.W  # h=6km

    # maximum angle of climb: 1) solve for the velocity for maximum gamma
    poly_s = [
        1 / 2 * cp1data.rho_s ** 2 * cp1data.C_D_0 / cp1data.WS,
        0,
        0,
        1 / 2 * cp1data.rho_s * cp1data.PAs / cp1data.W
        - 2 * cp1data.K * cp1data.WS,
    ]
    poly_h = [
        1 / 2 * rho_h ** 2 * cp1data.C_D_0 / cp1data.WS,
        0,
        0,
        1 / 2 * rho_h * (cp1data.PAs * sqrt(rho_h / cp1data.rho_s)) / cp1data.W
        - 2 * cp1data.K * cp1data.WS,
    ]
    Vgamma_s = np.roots(poly_s)
    Vgamma_h = np.roots(poly_h)
    Vgammamax_s_est = (
        2
        * cp1data.K
        * cp1data.WS
        / (1 / 2 * cp1data.rho_s * cp1data.PAs / cp1data.W)
    )
    Vgammamax_h_est = (
        2
        * cp1data.K
        * cp1data.WS
        / (
            1
            / 2
            * rho_h
            * (cp1data.PAs * sqrt(rho_h / cp1data.rho_s))
            / cp1data.W
        )
    )

    # maximum angle of climb: 2) compute gamma from the velocity (sea level)
    q_gamma_s = 1 / 2 * cp1data.rho_s * Vgammamax_s_est ** 2
    CL_gamma_s = cp1data.W / (q_gamma_s * cp1data.S)
    CD_gamma_s = cp1data.C_D_0 + cp1data.K * CL_gamma_s ** 2
    D_gamma_s = q_gamma_s * cp1data.S * CD_gamma_s
    RC_gamma_s = (cp1data.PAs - D_gamma_s * Vgammamax_s_est) / cp1data.W
    gammamax_s = np.arcsin(RC_gamma_s / Vgammamax_s_est) * 180 / pi

    # maximum angle of climb: 2) compute gamma from the velocity (h=6km)
    q_gamma_h = 1 / 2 * cp1data.rho_s * Vgammamax_h_est ** 2
    CL_gamma_h = cp1data.W / (q_gamma_h * cp1data.S)
    CD_gamma_h = cp1data.C_D_0 + cp1data.K * CL_gamma_h ** 2
    D_gamma_h = q_gamma_h * cp1data.S * CD_gamma_h
    RC_gamma_h = (
        cp1data.PAs * sqrt(rho_h / cp1data.rho_s) - D_gamma_h * Vgammamax_h_est
    ) / cp1data.W
    gammamax_h = np.arcsin(RC_gamma_h / Vgammamax_h_est) * 180 / pi

    ## power available to sustain given speeds
    P = cp1data.PAs * sqrt(rho_h / cp1data.rho_s)
    v = np.arange(0.1, 100, 1)
    p1 = [
        0.5 * rho_h * cp1data.S * cp1data.C_D_0,
        0,
        0,
        -P,
        cp1data.K * cp1data.S * cp1data.WS ** 2 / (0.5 * rho_h),
    ]
    f1 = np.polyval(p1, v)

    plt.figure(2)
    plt.plot(v, f1, linewidth=2)
    plt.grid()
    plt.title("speed profile under P_A = {:.4}".format(P))
    plt.xlabel("velocity (m/cp1data.s)")
    plt.ylabel("P_A (N m/cp1data.s)")

    plt.show()


if __name__ == "__main__":
    main()