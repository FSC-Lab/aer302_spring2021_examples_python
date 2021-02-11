import SAM
from aircraft import CP1
from funcs import *

import numpy as np
import matplotlib.pyplot as plt

h = 6000

# performance at sea level
Vminpr = 60
PR_sea = power_required(CP1, Vminpr, SAM.air_sea)
PA_sea = power_available(CP1, SAM.air_sea)
RC_sea, gmm_sea = climb_rate(CP1, Vminpr, SAM.air_sea)

Vs = np.linspace(10, 210, 100)
PR_Vs, PA_Vs, RC_Vs, Gmm_Vs = power_vs_v(
								CP1, SAM.air_sea, Vs)
Vsx, Vsy = hodograph(Vs, Gmm_Vs)


# performance at h
Vh = np.linspace(10, 210, 100)
PR_Vh, PA_Vh, RC_Vh, Gmm_Vh = power_vs_v(
								CP1, SAM.tropo(h), Vh)
V_PRmin_h = min_power_required(CP1, SAM.tropo(h))
RCmax_h = climb_rate(CP1, V_PRmin_h, SAM.tropo(h))
Vhx, Vhy = hodograph(Vh, Gmm_Vh)


# TODO: polish plots
plt.plot(Vs, PA_Vs)
plt.plot(Vs, PR_Vs)
plt.plot(Vh, PA_Vh)
plt.plot(Vh, PR_Vh)
plt.grid()
# plt.plot(Vsx, Vsy)
plt.show()