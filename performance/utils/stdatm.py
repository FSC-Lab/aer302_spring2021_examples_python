# stdatm.m - standard atmosphere (SI units)
# def [Temperature (K), pressure (N/m**2), density (kg/m**3)] = stdatm(altitude [m])
#
# created on: 25-Sep-00
# updated on: 27-Sep-00
from math import exp
from typing import Tuple

def stdatm(h: float) -> Tuple[float, float, float]:
    """Standard Atmosphere Model

    Args:
         h (float): Altitude in meters, must be positive

    Returns:
         Tuple[float, float, float]: A tuple of ambient temperature, standard pressure, and air pressure at the given altitude
    """

    # parameters
    a = -6.5e-03  # gradiant slope
    g = 9.8
    R = 287

    # sea level constants
    rho_s = 1.2250  # kg/m**3
    p_s = 1.01325e05  # N /m**2
    T_s = 288.16  # K

    # 11km constants
    T_11 = T_s + a * 11000
    p_11 = p_s * (T_11 / T_s) ** (-g / (a * R))
    rho_11 = rho_s * (T_11 / T_s) ** (-(g / (a * R) + 1))

    # temperature on altitude h
    if h < 0:
        print("Expected positive altitude value, got".format(h))
    elif h <= 11000:
        T = T_s + a * h
        p = p_s * (T / T_s) ** (-g / (a * R))
        rho = rho_s * (T / T_s) ** (-(g / (a * R) + 1))
    elif h > 11000:
        T = T_11
        p = p_11 * exp(-g / (R * T) * (h - 11000))
        rho = rho_11 * exp(-g / (R * T) * (h - 11000))

    return T, p, rho
