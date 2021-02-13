#!/usr/bin/python3

# cj1ac.m - absolute ceiling and time to climb
#
# created on: 25-Sep-00
# updated on: 27-Sep-00
#

from math import sqrt

import matplotlib.pyplot as plt
import numpy as np

# run data file first
from data import cj1data
from utils import stdatm


def main():
    # altitude 0 - 25000 m
    h = np.arange(0, 25000, 50)

    T, p, rho = stdatm(h)

    sigma: np.ndarray = rho / cj1data.rho_s  # density ratio
    # calculate thrust load
    TW = cj1data.F_s * sigma / cj1data.W
    # V for RC max
    Vrcmax = np.sqrt(
        cj1data.WS
        * TW
        / (3 * rho * cj1data.C_D_0)
        * (1 + np.sqrt(1 + 12 * cj1data.K * cj1data.C_D_0 / TW ** 2))
    )
    # max rate of climb
    q_rcmax = 1 / 2 * rho * Vrcmax ** 2

    RCmax_h = Vrcmax * (
        TW
        - q_rcmax * cj1data.C_D_0 / cj1data.WS
        - cj1data.K * cj1data.WS / q_rcmax
    )

    plt.figure(1)
    plt.plot(RCmax_h, h, "-", linewidth=2)
    plt.grid()
    plt.axis([0, 45, 0, 20000])
    plt.title("ceiling")
    plt.xlabel(" RC max (m/s)")
    plt.ylabel(" altitude (m)")
    plt.show()
    # plotdlg

    # plot(h,InvMaxRC,'-')
    # title('time to climb')
    # ylabel(' 1 / RC max (m/s)')
    # xlabel(' altitude (m)')
    # axis([0 6100 0 0.1])
    # plotdlg
    #


if __name__ == "__main__":
    main()
