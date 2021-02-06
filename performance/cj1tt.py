# cj1tt.m - thrust  at sea level
#
# created on: 24-Sep-00
# updated on:
#

## basic data for CP1, P45 lecture 08
from math import pi, sqrt

import matplotlib.pyplot as plt
import numpy as np

## basic data for CP1, P16 lecture 08
import data
from utils import stdatm

def main():
    # some data at altitude h=6.7km
    h = 6705.6
    T_h, p_h, rho_h = stdatm(h)

    ## compute speed for given thrust available (P47 lecture 08)
    # at sea level
    c_s = [
        1 - (data.F_s / data.S) / data.C_D_0,
        data.K * data.WS ** 2 / data.C_D_0,
    ]
    q_s = np.roots(c_s)
    V_s = sqrt(2 * q_s / data.rho_s)
    Vstall_s = sqrt(2 * data.W / (data.rho_s * data.S * data.C_L_max))

    # at altitude h=6.7km
    F_h = (
        rho_h / data.rho_s * data.F_s
    )  # thrust available at h is different from sea  level
    c_h = [1 - (F_h / data.S) / data.C_D_0, data.K * data.WS ** 2 / data.C_D_0]
    q_h = np.roots(c_h)
    V_h = sqrt(2 * q_h / rho_h)
    Vstall_h = sqrt(
        2 * data.W / (rho_h * data.S * data.C_L_max)
    )  # assume C_L_max indepent of the altitude

    ## thrut required and thrust available plt.plots

    # at sea level
    Vs = np.arange(10, 410)
    qsS = 0.5 * data.rho_s * Vs ** 2 * data.S
    C_Ls = data.W / qsS
    C_Ds = data.C_D_0 + data.K * C_Ls ** 2
    Ds = qsS * C_Ds  # thrust required
    Fs = data.F_s  # thrust available

    # at altitude h=6.7km
    Vh = Vs * sqrt(data.rho_s / rho_h)
    qS = 0.5 * rho_h * Vh ** 2 * data.S
    C_L = data.W / qS
    C_D = data.C_D_0 + data.K * C_L ** 2
    Dh = qS * C_D  # thrust required
    Fh = (rho_h / data.rho_s) * Fs  # thrust available


    # minimum Thrust
    TRmin_s, iMin_s = (min(Ds), np.argmin(Ds))
    VTRmin_s = 10 + iMin_s * 1
    TRmin_h, iMin_h = (min(Dh), np.argmin(Dh))
    VTRmin_h = (10 + iMin_h * 1) * sqrt(data.rho_s / rho_h)

    # plt.plot
    plt.figure(1)
    plt.plot(Vs, Ds, "-", label="T_R: sea-level")
    plt.axhline(y=Fs, linestyle="-", label="T_A: sea-level")
    plt.plot(Vh, Dh, "--", label="T_R: altitude")
    plt.axhline(y=Fh, linestyle="--", label="T_A: altitude")
    plt.plot(VTRmin_s, TRmin_s, "o", VTRmin_h, TRmin_h, "*")
    plt.grid()
    plt.title("Thrust Required and Available")
    plt.xlabel("velocity (m/s)")
    plt.ylabel("thrust (N)")

    ## CL/CD progile, illustration @ sea level
    LDs = C_Ls / C_Ds
    LDmax_s = max(LDs)
    CLmax_s = sqrt(data.C_D_0 / data.K)

    # lift to drag ratio vs. drag
    plt.figure(2)
    plt.plot(C_Ls, LDs, CLmax_s, LDmax_s, "o")
    plt.axis([0, 5, 0, 20])
    plt.grid()
    plt.title("Lift to Drag Ratio")
    plt.xlabel("C_L")
    plt.ylabel("C_L / C_D")

    # drag pola, lift coefficient vs. drag coefficient
    plt.figure(3)
    plt.plot(C_L, C_D)
    plt.axis([0, 5, 0, 1.5])
    plt.grid()
    plt.title("Drag Polar")
    plt.xlabel("C_L")
    plt.ylabel("C_D")

    plt.show()

if __name__ == "__main__":
    main()