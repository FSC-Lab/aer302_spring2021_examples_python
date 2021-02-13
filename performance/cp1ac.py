#!/usr/bin/env python3
# cp1ac.m - absolute and service ceilings
#
# created on: 10.Feb.2020
# updated on:
#

# run data file first
from math import pi

import matplotlib.pyplot as plt
import numpy as np

from data import cp1data
from utils import safe_asin, stdatm


def main():
    h = np.arange(0, 25000, 50)  # altitude 0 - 25000 m
    T, p, rho = stdatm(h)

    sigma = np.sqrt(rho / cp1data.rho_s)  # density ratio
    # calculate power
    PA_h = cp1data.PAs * sigma
    # min PR
    PRmin = (
        cp1data.W
        * np.sqrt(cp1data.K * cp1data.C_D_0)
        * 4
        / np.sqrt(3)
        * np.sqrt(
            2 * cp1data.WS / np.sqrt(3 * cp1data.C_D_0 / cp1data.K) / rho
        )
    )

    # max RC
    RCmax_h = PA_h / cp1data.W - PRmin / cp1data.W

    plt.figure(1)
    plt.plot(RCmax_h, h, "-", linewidth=2)
    plt.axis([0, 10, 0, 15000])
    plt.title("ceiling")
    plt.xlabel(" RC max (m/s)")
    plt.ylabel(" altitude (m)")

    plt.show()


# plt.plot(h,InvMaxRC,'-')
# plt.title('time to climb')
# plt.ylabel(' 1 / RC max (m/s)')
# plt.xlabel(' altitude (m)')
# plt.axis([0 6100 0 0.1])
# plt.plotdlg

if __name__ == "__main__":
    main()