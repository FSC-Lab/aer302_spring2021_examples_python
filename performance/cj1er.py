#!/usr/bin/env python3
#
# cj1er.m - endurance and range
#
# created on: 24-Sep-00
# updated on:
#

import argparse

# run data file first
from math import pi, sqrt

import matplotlib.pyplot as plt
import numpy as np

from data import cj1data
from utils import safe_asin, stdatm


def main():
    ap = argparse.ArgumentParser("endurance and range")
    ap.add_argument(
        "-v", "--verbose", action="store_true", help="Toggles verbose output"
    )

    args = ap.parse_args()
    # altitude
    h = 6705.6  # (m)
    T, p, rho = stdatm(h)

    # lift-to-drag ratio
    V = np.arange(10, 410)
    # at sea level
    C_L_s = cj1data.W / (0.5 * cj1data.rho_s * V ** 2 * cj1data.S)
    C_D_s = cj1data.C_D_0 + cj1data.K * C_L_s ** 2
    LD_s = C_L_s / C_D_s
    LsD_s = np.sqrt(C_L_s) / C_D_s
    # at altitude
    C_L = cj1data.W / (0.5 * rho * V ** 2 * cj1data.S)
    C_D = cj1data.C_D_0 + cj1data.K * C_L ** 2
    LD = C_L / C_D
    LsD = np.sqrt(C_L) / C_D

    # maximum ratio value
    LD_sMax = np.max(LD_s)
    LsD_sMax = np.max(LsD_s)
    LDMax = np.max(LD)
    LsDMax = np.max(LsD)

    # calculation of endurance and range
    W0 = cj1data.W
    W1 = cj1data.W - cj1data.Wf  # empty tank
    Es = 1 / cj1data.c_t * LD_sMax * np.log(W0 / W1) / 3600  # hour
    E = 1 / cj1data.c_t * LDMax * np.log(W0 / W1) / 3600  # hour
    Rs = (
        2
        * np.sqrt(2 / (rho * cj1data.S))
        * 1
        / cj1data.c_t
        * LsD_sMax
        * (np.sqrt(W0) - np.sqrt(W1))
    )
    R = (
        2
        * np.sqrt(2 / (rho * cj1data.S))
        * 1
        / cj1data.c_t
        * LsDMax
        * (np.sqrt(W0) - np.sqrt(W1))
    )

    # presentation of results
    plt.plot(V, LD_s, "-", label="C_L/C_D sea level", linewidth=2)
    plt.plot(V, LsD_s, "-", label="\sqrt(C_L)/C_D sea level", linewidth=2)
    plt.plot(V, LD, "--", label="C_L/C_D h", linewidth=2)
    plt.plot(V, LsD, "--", label="\sqrt(C_L)/C_D h", linewidth=2)
    plt.grid
    plt.title("C_L to C_D ratio")
    plt.xlabel(" velocity (m/cj1data.s)")
    plt.ylabel(" ")

    plt.show()

    if args.verbose:

        msg = (
            "- the maximum C_L/C_D calcualted from sea level is:     {:^12.4}"
            "- the endurance is:                                     {:^12.4}(hour)"
            "- the maximum C_L**1/2/C_D calcualted from sea level is:{:^12.4}"
            "- the range is:                                         {:^12.4}(m)"
            "- the maximum C_L/C_D calcualted from altitude  is:     {:^12.4}"
            "- the endurance is:                                     {:^12.4}(hour)"
            "- the maximum C_L**1/2/C_D calcualted from altitude  is:{:^12.4}"
            "- the range is:                                         {:^12.4}(m)"
        )

        print(msg.format(LD_sMax, Es, LsD_sMax, Rs, LDMax, E, LsDMax, R))


if __name__ == "__main__":
    main()
