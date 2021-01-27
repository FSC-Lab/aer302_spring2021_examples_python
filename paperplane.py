import argparse
from math import atan, cos, pi, sin, sqrt

import matplotlib.pyplot as plt
import numpy as np
from control import initial_response, ss, step_response
from scipy.integrate import ode

ap = argparse.ArgumentParser("AER302 Paperplane example")
ap.add_argument(
    "--latex",
    action="store_true",
    default=False,
    help="Toggles use of LaTeX backend on / off",
)

args = ap.parse_args()
if args.latex:
    from matplotlib import rc

    rc("text", usetex=True)

""" basic parameters """
S = 0.017  # Reference Area, m^2
AR = 0.86  # Wing Aspect Ratio
e = 0.9  # Oswald Efficiency Factor
m = 0.003  # Mass, kg
g = 9.8  # Gravitational acceleration, m/s^2
rho = 1.225  # Air density at Sea Level, kg/m^3
CLa = (
    pi * AR / (1 + sqrt(1 + (AR / 2) ** 2))
)  # Lift-Coefficient Slope, per rad
CDo = 0.02  # Zero-Lift Drag Coefficient
kappa = 1 / (3.141592 * e * AR)  # Induced Drag Factor

""" optimal Condition at (CL/CD)max """
CL_star = sqrt(CDo / kappa)  # CL for Maximum Lift/Drag Ratio
CD_star = CDo + kappa * CL_star ** 2  # Corresponding CD
LD_star = CL_star / CD_star  # Maximum Lift/Drag Ratio
Gamma_star = -atan(1 / LD_star)  # Corresponding Flight Path Angle, rad
Gamma_starDeg = Gamma_star * 180 / pi
V_star = sqrt(2 * m * g * cos(Gamma_star) / (rho * S * (CL_star)))
# Ve = sqrt(2 * m * g /(rho * S * (CLe * cos(Gammae) - CDe * sin(Gammae))))
# Corresponding Velocity, m/s
Alpha_star = CL_star / CLa  # Corresponding Angle of Attack, rad
Alpha_starDeg = Alpha_star * 180 / pi

""" equilibrium condition  under a fixed alpha """
Alphae = 8 * pi / 180  # assume constant angle of attack, rad
# Alpha =   Alphae + 10 * pi/180		# angle of attack, rad
CLe = CLa * Alphae
CDe = CDo + kappa * CLe ** 2
CL = CLe
CD = CDe
Gammae = -atan(1 / (CLe / CDe))
Ve = sqrt(2 * m * g * cos(Gammae) / (rho * S * (CLe)))

""" initial conditions """
Vo = Ve  # initial velocity, m/s
H = 2  # Initial Height, m
R = 0  # Initial Range, m
Gammao = Gammae  # Initial flight path angle, rad
to = 0  # Initial Time, sec
tf = 6  # Final Time, sec
tspan = [to, tf]

""" nolinear system solution """


def EqMotion(t, x):
    # Fourth-Order Equations of Aircraft Motion

    V = x[0]
    Gam = x[1]
    q = 0.5 * rho * V ** 2  # Dynamic Pressure, N/m^2

    xdot = np.array(
        [
            (-CD * q * S - m * g * sin(Gam)) / m,
            (CL * q * S - m * g * cos(Gam)) / (m * V),
            V * sin(Gam),
            V * cos(Gam),
        ]
    )

    return xdot


# integrator


def integrate(xo, to=to, tf=tf, dt=0.1):

    solver = ode(EqMotion).set_integrator("dopri5")
    # a) Equilibrium Glide at Maximum Lift/Drag Ratio
    ts, xs = [], []
    solver.set_initial_value(xo, to)
    while solver.successful() and solver.t < tf:
        solver.integrate(solver.t + 0.1)
        ts.append(solver.t)
        xs.append(solver.y)

    return np.asarray(ts), np.asarray(xs)


# a) Equilibrium Glide at Maximum Lift/Drag Ratio
ta, xa = integrate(np.array([Ve, Gammae, H, R]))

# # b) Oscillating Glide due to Zero Initial Flight Path Angle
tb, xb = integrate(np.array([Ve, 0, H, R]))

# c) Effect of Increased Initial Velocity
tc, xc = integrate(np.array([1.5 * Ve, 0, H, R]))

# d) Effect of Further Increase in Initial Velocity
td, xd = integrate(np.array([3.0 * Ve, 0, H, R]))


""" state-space model and solution """
# no control input
A = np.array(
    [
        [-rho * Ve * CD * S / m, -g * cos(Gammae)],
        [rho * CL * S / m, g * sin(Gammae) / Ve],
    ]
)
B = np.zeros((2, 1))
C = np.identity(2)
D = np.zeros((2, 1))
sys1 = ss(A, B, C, D)

X0a = np.array([Ve, Gammae])
X0b = np.array([Ve, 0])
X0c = np.array([1.5 * Ve, 0])
X0d = np.array([3 * Ve, 0])

tl = np.linspace(to, tf, int((tf - to) / 0.01))
_, ya = initial_response(sys1, X0=X0a, T=tl)
_, yb = initial_response(sys1, X0=X0b, T=tl)
_, yc = initial_response(sys1, X0=X0c, T=tl)
_, yd = initial_response(sys1, X0=X0d, T=tl)


# step control input
Ahat = np.array(
    [
        [2 * g * sin(Gammae) / Ve, -g * cos(Gammae)],
        [2 * g * cos(Gammae) / (Ve ** 2), g * sin(Gammae) / Ve],
    ]
)
Bhat = np.array(
    [[2 * g * cos(Gammae) * kappa * CLa], [g * cos(Gammae) / (Ve * Alphae)]]
)
Chat = np.identity(2)
Dhat = np.zeros((2, 1))
sys2 = ss(Ahat, Bhat, Chat, Dhat)

deltaAlpha = 10 * pi / 180
_, y = step_response(sys2, T=tl)
Vl = y[0, :] + Ve
Gammal = y[1, :] * deltaAlpha + Gammae

""" plot """
lw = 2
fs = 20

# Longitudinal trajectory
plt.title("Longitudinal Trajectory", fontsize=fs)

plt.plot(xa[:, 3], xa[:, 0], lw=lw, label=r"$V_e,\gamma_e$")
plt.plot(xb[:, 3], xb[:, 0], lw=lw, label=r"$V_e, \gamma_0=0$")
plt.plot(xc[:, 3], xc[:, 0], lw=lw, label=r"$1.5V_e, \gamma_0=0$")
plt.plot(xd[:, 3], xd[:, 0], lw=lw, label=r"$3V_e, \gamma_0=0$")

plt.legend(fontsize=fs)
plt.xlabel("Distance, " + r"$m$", fontsize=fs)
plt.ylabel("Height, " + r"$m$", fontsize=fs)
plt.grid()
plt.gca().tick_params(axis="both", which="major", labelsize=fs)
plt.gca().tick_params(axis="both", which="minor", labelsize=fs)
# plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.15)
plt.tight_layout()

plt.show()


# State variables
fig, axs = plt.subplots(2, 2, figsize=(10, 10))

axs[0, 0].plot(ta, xa[:, 0], lw=lw, label=r"$V_e,\gamma_e$")
axs[0, 0].plot(tb, xb[:, 0], lw=lw, label=r"$V_e, \gamma_0=0$")
axs[0, 0].plot(tc, xc[:, 0], lw=lw, label=r"$1.5V_e, \gamma_0=0$")
axs[0, 0].plot(td, xd[:, 0], lw=lw, label=r"$3V_e, \gamma_0=0$")
axs[0, 0].set_xlabel("Time, " + r"$t$", fontsize=fs)
axs[0, 0].set_ylabel("Velocity, " + r"$m/s$", fontsize=fs)

axs[0, 1].plot(ta, xa[:, 1], lw=lw, label=r"$V_e,\gamma_e$")
axs[0, 1].plot(tb, xb[:, 1], lw=lw, label=r"$V_e, \gamma_0=0$")
axs[0, 1].plot(tc, xc[:, 1], lw=lw, label=r"$1.5V_e, \gamma_0=0$")
axs[0, 1].plot(td, xd[:, 1], lw=lw, label=r"$3V_e, \gamma_0=0$")
axs[0, 1].set_xlabel("Time, " + r"$t$", fontsize=fs)
axs[0, 1].set_ylabel("Flight Path Angle, " + r"$rad$", fontsize=fs)

axs[1, 0].plot(ta, xa[:, 2], lw=lw, label=r"$V_e,\gamma_e$")
axs[1, 0].plot(tb, xb[:, 2], lw=lw, label=r"$V_e, \gamma_0=0$")
axs[1, 0].plot(tc, xc[:, 2], lw=lw, label=r"$1.5V_e, \gamma_0=0$")
axs[1, 0].plot(td, xd[:, 2], lw=lw, label=r"$3V_e, \gamma_0=0$")
axs[1, 0].set_xlabel("Time, " + r"$t$", fontsize=fs)
axs[1, 0].set_ylabel("Altitude, " + r"$m$", fontsize=fs)

axs[1, 1].plot(ta, xa[:, 3], lw=lw, label=r"$V_e,\gamma_e$")
axs[1, 1].plot(tb, xb[:, 3], lw=lw, label=r"$V_e, \gamma_0=0$")
axs[1, 1].plot(tc, xc[:, 3], lw=lw, label=r"$1.5V_e, \gamma_0=0$")
axs[1, 1].plot(td, xd[:, 3], lw=lw, label=r"$3V_e, \gamma_0=0$")
axs[1, 1].set_xlabel("Time, " + r"$t$", fontsize=fs)
axs[1, 1].set_ylabel("Range, " + r"$rad$", fontsize=fs)

for ax in axs:
    for axx in ax:
        axx.grid()
        axx.tick_params(axis="both", which="major", labelsize=fs)
        axx.tick_params(axis="both", which="minor", labelsize=fs)
axs[0, 1].legend(fontsize=fs)
plt.tight_layout()

plt.show()

# no control vs. control based on linearized system
fig, axs = plt.subplots(1, 2, figsize=(10, 5))

axs[0].plot(ta, xa[:, 0], lw=lw, label="no control")
axs[0].plot(
    tl, Vl, lw=lw, label=r"linearized control,\n$\Delta\alpha = 10 deg$"
)
axs[0].set_xlabel("Time, " + r"$t$", fontsize=fs)
axs[0].set_ylabel("Velocity, " + r"$m/s$", fontsize=fs)

axs[1].plot(ta, xa[:, 0], lw=lw, label="no control")
axs[1].plot(
    tl, Gammal, lw=lw, label=r"linearized control,\n$\Delta\alpha = 10 deg$"
)
axs[1].set_xlabel("Time, " + r"$t$", fontsize=fs)
axs[1].set_ylabel("Flight Path Angle, " + r"$rad$", fontsize=fs)

for ax in axs:
    ax.grid()
    ax.tick_params(axis="both", which="major", labelsize=fs)
    ax.tick_params(axis="both", which="minor", labelsize=fs)
axs[1].legend(fontsize=fs)
plt.tight_layout()

plt.show()

# nonlinear vs. linear initial responses: velocity
fig, axs = plt.subplots(2, 2, figsize=(10, 10))

axs[0, 0].plot(ta, xa[:, 0], lw=lw, label="numerical")
axs[0, 0].plot(tl, ya[0, :], lw=lw, label="linearized")
axs[0, 0].set_xlabel("Time, " + r"$t$", fontsize=fs)
axs[0, 0].set_ylabel("Velocity, " + r"$m/s$", fontsize=fs)

axs[0, 1].plot(tb, xb[:, 0], lw=lw, label="numerical")
axs[0, 1].plot(tl, yb[0, :], lw=lw, label="linearized")
axs[0, 1].set_xlabel("Time, " + r"$t$", fontsize=fs)
axs[0, 1].set_ylabel("Velocity, " + r"$m/s$", fontsize=fs)

axs[1, 0].plot(tc, xc[:, 0], lw=lw, label="numerical")
axs[1, 0].plot(tl, yc[0, :], lw=lw, label="linearized")
axs[1, 0].set_xlabel("Time, " + r"$t$", fontsize=fs)
axs[1, 0].set_ylabel("Velocity, " + r"$m/s$", fontsize=fs)

axs[1, 1].plot(td, xd[:, 0], lw=lw, label="numerical")
axs[1, 1].plot(tl, yd[0, :], lw=lw, label="linearized")
axs[1, 1].set_xlabel("Time, " + r"$t$", fontsize=fs)
axs[1, 1].set_ylabel("Velocity, " + r"$m/s$", fontsize=fs)

for ax in axs:
    for axx in ax:
        axx.grid()
        axx.tick_params(axis="both", which="major", labelsize=fs)
        axx.tick_params(axis="both", which="minor", labelsize=fs)
axs[0, 1].legend(fontsize=fs)
plt.tight_layout()

plt.show()

# nonlinear vs. linear initial responses: flight path angle
fig, axs = plt.subplots(2, 2, figsize=(10, 10))

axs[0, 0].plot(ta, xa[:, 1], lw=lw, label="numerical")
axs[0, 0].plot(tl, ya[1, :], lw=lw, label="linearized")
axs[0, 0].set_xlabel("Time, " + r"$t$", fontsize=fs)
axs[0, 0].set_ylabel(r"$\gamma$", fontsize=fs)

axs[0, 1].plot(tb, xb[:, 1], lw=lw, label="numerical")
axs[0, 1].plot(tl, yb[1, :], lw=lw, label="linearized")
axs[0, 1].set_xlabel("Time, " + r"$t$", fontsize=fs)
axs[0, 1].set_ylabel(r"$\gamma$", fontsize=fs)

axs[1, 0].plot(tc, xc[:, 1], lw=lw, label="numerical")
axs[1, 0].plot(tl, yc[1, :], lw=lw, label="linearized")
axs[1, 0].set_xlabel("Time, " + r"$t$", fontsize=fs)
axs[1, 0].set_ylabel(r"$\gamma$", fontsize=fs)

axs[1, 1].plot(td, xd[:, 1], lw=lw, label="numerical")
axs[1, 1].plot(tl, yd[1, :], lw=lw, label="linearized")
axs[1, 1].set_xlabel("Time, " + r"$t$", fontsize=fs)
axs[1, 1].set_ylabel(r"$\gamma$", fontsize=fs)

for ax in axs:
    for axx in ax:
        axx.grid()
        axx.tick_params(axis="both", which="major", labelsize=fs)
        axx.tick_params(axis="both", which="minor", labelsize=fs)
axs[0, 1].legend(fontsize=fs)
plt.tight_layout()

plt.show()
