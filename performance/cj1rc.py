#!/usr/bin/env python3

# cj1rc.m - rate of climbing
#
# created on: 03-Oct-00
# updated on:
#
# Revision History


# run data file first
from math import pi, sqrt

import matplotlib.pyplot as plt
import numpy as np

## basic data for CP1, P16 lecture 08
from data import cj1data
from utils import stdatm, safe_asin


def main():
    # altitude
    h = 6705.6
    T_h, p_h, rho_h = stdatm(h)
    F_h = (rho_h / cj1data.rho_s) * cj1data.F_s  # TA

    # RC at sea level - assume gamma (flight path angle) << 1
    Vs = np.arange(10, 400)

    # at sea level
    qsS = 0.5 * cj1data.rho_s * Vs ** 2 * cj1data.S
    C_Ls = cj1data.W / qsS
    C_Ds = cj1data.C_D_0 + cj1data.K * C_Ls ** 2
    Ds = qsS * C_Ds  # TR
    Fs = cj1data.F_s  # TA
    RC_s = (Fs * Vs - Ds * Vs) / cj1data.W
    gamma_s = safe_asin(RC_s / Vs, deg=True)
    # at altitude
    Vh = Vs * sqrt(cj1data.rho_s / rho_h)
    qS = 0.5 * rho_h * Vh ** 2 * cj1data.S
    C_L = cj1data.W / qS
    C_D = cj1data.C_D_0 + cj1data.K * C_L ** 2
    Dh = qS * C_D  # TR

    Fh = F_h  # TA
    RC_h = (Fh * Vh - Dh * Vh) / cj1data.W
    gamma_h = safe_asin(RC_h / Vh, deg=True)

    # compared to power curve
    PR_s = Ds * Vs
    PA_s = cj1data.F_s * Vs

    # hodograph (approximation)
    Vv_s = RC_s
    # Vh_s = sqrt(V**2 - Vv_s**2)
    Vh_s = Vs * np.cos(np.deg2rad(gamma_s))
    Vv_h = RC_h
    Vh_h = Vh * np.cos(np.deg2rad(gamma_h))

    # compare maximum values
    RCmax_s = np.max(RC_s)
    Vrcmax_s = Vs[RC_s == RCmax_s]
    gammamax_s = np.max(gamma_s)
    Vgammamax_s = Vs[gamma_s == gammamax_s]

    # analytical calculation of maximum values P10 lecture 09
    TWs = cj1data.F_s / cj1data.W
    Vrcmax_s_calc = sqrt(
        (cj1data.WS * TWs)
        / (3 * cj1data.rho_s * cj1data.C_D_0)
        * (1 + sqrt(1 + 12 * cj1data.K * cj1data.C_D_0 / TWs ** 2))
    )

    # alternative calculation use V* as base P11 lecture 09
    Vstar = sqrt(
        2 * cj1data.WS / (cj1data.rho_s * sqrt(cj1data.C_D_0 / cj1data.K))
    )
    z = sqrt(TWs ** 2 / (12 * cj1data.K * cj1data.C_D_0) + 1)
    Vrcmax_s_calc2 = (
        Vstar * sqrt(1 / (2 * sqrt(3))) * (sqrt(z + 1) + sqrt(z - 1))
    )
    Vgammamax_s_calc = sqrt(
        2 * cj1data.WS / (cj1data.rho_s * sqrt(cj1data.C_D_0 / cj1data.K))
    )

    # Eqn = [-cj1data.rho_s**2 * cj1data.S**2 * cj1data.C_d_0, cj1data.rho_s * cj1data.S * cj1data.F_s, 0, 2*cj1data.K*cj1data.W**2]
    # Vcal = roots(Eqn)

    fig_1, axs = plt.subplots(2, 1)

    axs[0].plot(Vs, PA_s, "-", label="PA", linewidth=2)
    axs[0].plot(Vs, PR_s, "--", label="PR", linewidth=2)
    axs[0].grid()
    axs[0].set_title("Power available and required at Sea-Level")
    axs[0].set_xlabel("velocity (m/cj1data.s)")
    axs[0].set_ylabel("Power (N.m/cj1data.s)")
    axs[0].legend()

    axs[1].plot(Vs, RC_s, "-", linewidth=2)
    axs[1].plot(Vrcmax_s, RCmax_s, "o", linewidth=2)
    axs[1].grid()
    axs[1].axis([0, 400, -100, 50])
    axs[1].set_title("Rate of Climb at Sea-Level")
    axs[1].set_xlabel("velocity (m/cj1data.s)")
    axs[1].set_ylabel("R/C (m/cj1data.s)")

    fig_2, axs = plt.subplots(2, 1)
    axs[0].plot(Vs, PA_s, "-", label="PA", linewidth=2)
    axs[0].plot(Vs, PR_s, "--", label="PR", linewidth=2)
    axs[0].grid()
    axs[0].set_title("Power available and required at Sea-Level")
    axs[0].set_xlabel("velocity (m/cj1data.s)")
    axs[0].set_ylabel("Power (N.m/cj1data.s)")

    axs[1].plot(Vs, gamma_s, "-", linewidth=2)
    axs[1].plot(Vgammamax_s, gammamax_s, "o", linewidth=2)
    axs[1].grid()
    axs[1].axis([0, 400, -100, 50])
    axs[1].set_title("Climbing Angle at Sea-Level")
    axs[1].set_xlabel("velocity (m/cj1data.s)")
    axs[1].set_ylabel("\gamma (degree)")

    plt.figure(3)
    plt.plot(Vh_s, Vv_s, "-", label="T_A ={}".format(cj1data.F_s))
    plt.plot(Vh_h, Vv_h, "--", label="T_A={}".format(F_h))
    plt.grid()
    plt.title("Hodograph for [CJ1]")
    plt.axis([0, 400, 0, 50])
    plt.xlabel("Vh (m/cj1data.s)")
    plt.ylabel("Vv (m/cj1data.s)")
    plt.legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()